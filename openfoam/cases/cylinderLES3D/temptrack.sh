#!/usr/bin/env bash
set -euo pipefail

INTERVAL="${1:-2}"              # Sekunden zwischen Messungen
LOG_DIR="${XDG_STATE_HOME:-$HOME/.local/state}/temptrack"
LOG_FILE="$LOG_DIR/temps.csv"

mkdir -p "$LOG_DIR"

# --- Helfer: CPU Temp (best effort) ---
get_cpu_temp() {
  # Versuch 1: "Package id 0" (Intel) oder ähnliches
  local t
  t="$(sensors 2>/dev/null | awk '
    /Package id 0:|Tctl:|Tdie:|CPU Temperature:/ {
      for (i=1;i<=NF;i++) if ($i ~ /^\+?[0-9]+(\.[0-9]+)?°C$/) { gsub(/[+°C]/,"",$i); print $i; exit }
    }')"
  if [[ -n "${t:-}" ]]; then
    echo "$t"; return
  fi

  # Versuch 2: irgendeinen "Core 0" nehmen
  t="$(sensors 2>/dev/null | awk '
    /Core 0:|Core0:/ {
      for (i=1;i<=NF;i++) if ($i ~ /^\+?[0-9]+(\.[0-9]+)?°C$/) { gsub(/[+°C]/,"",$i); print $i; exit }
    }')"
  if [[ -n "${t:-}" ]]; then
    echo "$t"; return
  fi

  echo "NaN"
}

# --- Helfer: NVIDIA GPU Temp (best effort) ---
get_gpu_temp() {
  if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=temperature.gpu --format=csv,noheader,nounits 2>/dev/null | head -n1 | tr -d ' '
    return
  fi
  echo "NaN"
}

# --- ASCII Balken ---
bar() {
  # bar <value> <max> <width>
  local v="$1" max="$2" width="$3"
  if [[ "$v" == "NaN" ]]; then
    printf "%-*s" "$width" "?"
    return
  fi
  # clamp
  awk -v v="$v" -v max="$max" -v w="$width" '
    BEGIN {
      if (v < 0) v=0
      if (v > max) v=max
      n = int((v/max)*w + 0.5)
      for (i=0; i<n; i++) printf "█"
      for (i=n; i<w; i++) printf " "
    }'
}

color_for_temp() {
  # Gibt ANSI-Farbe zurück je nach Temperatur
  local t="$1"
  [[ "$t" == "NaN" ]] && echo "0" && return
  awk -v t="$t" 'BEGIN{
    if (t < 60) print 32;      # grün
    else if (t < 80) print 33; # gelb
    else print 31;             # rot
  }'
}



# --- Init Log ---
if [[ ! -f "$LOG_FILE" ]]; then
  echo "timestamp_iso,cpu_c,gpu_c" > "$LOG_FILE"
fi

# --- Terminal vorbereiten ---
clear
tput civis
trap 'tput cnorm; echo; exit' INT TERM

echo "Logging to: $LOG_FILE"
echo "Interval: ${INTERVAL}s"
echo "Stop with: Ctrl+C"
echo

# Platzhalter-Zeilen (werden später überschrieben)
echo "CPU:"
echo "GPU:"
echo "Min/Max:"

# --- State ---
cpu_min=999; cpu_max=-1
gpu_min=999; gpu_max=-1


while true; do
  ts_iso="$(date -Iseconds)"
  ts="$(date +%H:%M:%S)"
  cpu="$(get_cpu_temp)"
  gpu="$(get_gpu_temp)"

  echo "${ts_iso},${cpu},${gpu}" >> "$LOG_FILE"

  # Min/Max (nur wenn Zahl)
  if [[ "$cpu" != "NaN" ]]; then
    cpu_i="${cpu%.*}"
    (( cpu_i < cpu_min )) && cpu_min="$cpu_i"
    (( cpu_i > cpu_max )) && cpu_max="$cpu_i"
  fi

  if [[ "$gpu" != "NaN" ]]; then
    gpu_i="${gpu%.*}"
    (( gpu_i < gpu_min )) && gpu_min="$gpu_i"
    (( gpu_i > gpu_max )) && gpu_max="$gpu_i"
  fi

  # Dynamische Balkenbreite je nach Terminalbreite
  cols="$(tput cols)"
  # grob: Text braucht ~38 Zeichen, Rest Balken
  bar_w=$(( cols - 45 ))
  (( bar_w < 10 )) && bar_w=10
  (( bar_w > 80 )) && bar_w=80

  cpu_bar="$(bar "$cpu" 100 "$bar_w")"
  gpu_bar="$(bar "$gpu" 100 "$bar_w")"

  cpu_col="$(color_for_temp "$cpu")"
  gpu_col="$(color_for_temp "$gpu")"


  # CPU Zeile (nach Header: wir sind bei Zeile 6 -> index 5)
  tput cup 5 0
  printf "CPU: %6s°C \033[%sm|%s|\033[0m  %s\033[K" "$cpu" "$cpu_col" "$cpu_bar" "$ts"

  # GPU Zeile
  tput cup 6 0
  printf "GPU: %6s°C \033[%sm|%s|\033[0m\033[K" "$gpu" "$gpu_col" "$gpu_bar"

  # Min/Max
  tput cup 7 0
  printf "Min/Max: CPU %s/%s °C   GPU %s/%s °C\033[K" \
    "${cpu_min/999/--}" "${cpu_max/-1/--}" "${gpu_min/999/--}" "${gpu_max/-1/--}"


  sleep "$INTERVAL"
done
