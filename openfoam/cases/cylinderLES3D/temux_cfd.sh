#!/usr/bin/env bash
set -euo pipefail

# ==========
# Settings
# ==========
SESSION="${1:-cfd}"
CASE_DIR="${2:-$PWD}"

# Commands (anpassen falls du andere Skripte nutzt)
SIM_CMD='echo "== Simulation =="; pwd; ls; echo; echo "Tipps:"; echo "  foamRun | mpirun -np 4 foamRun -parallel"; echo "  tail -f log.foamRun"; bash'
PARA_CMD='echo "== ParaView =="; pwd; echo; echo "Start: paraFoam"; bash'
TEMP_CMD='echo "== TempTrack =="; pwd; echo; echo "Start: ./temptrack.sh (oder tempreck.sh)"; bash'

# Optional: wenn du den Pfad/Name sicher kennst, z.B.
# TEMP_CMD='cd "'"$CASE_DIR"'"; ./temptrack.sh; bash'

# ==========
# tmux layout
# ==========
if tmux has-session -t "$SESSION" 2>/dev/null; then
  echo "Session '$SESSION' existiert schon. Attach mit: tmux attach -t $SESSION"
  exit 0
fi

tmux new-session -d -s "$SESSION" -c "$CASE_DIR" -n work

# links: Simulation
tmux send-keys -t "$SESSION:work.0" "cd '$CASE_DIR'; $SIM_CMD" C-m

# rechts: split vertical
tmux split-window -h -t "$SESSION:work" -c "$CASE_DIR"

# rechts unten: temp
tmux send-keys -t "$SESSION:work.1" "cd '$CASE_DIR'; $TEMP_CMD" C-m

# rechts oben: paraview (split rechts pane in zwei übereinander)
tmux split-window -v -t "$SESSION:work.1" -c "$CASE_DIR"
tmux send-keys -t "$SESSION:work.2" "cd '$CASE_DIR'; $PARA_CMD" C-m

# Größen hübsch machen: links ~60%, rechts ~40%
tmux resize-pane -t "$SESSION:work.0" -x 120 2>/dev/null || true

# Maus an (optional nice)
tmux set-option -t "$SESSION" -g mouse on

# Pane-Titel anzeigen
tmux set-option -t "$SESSION" -g pane-border-status top
tmux set-option -t "$SESSION" -g pane-border-format " #{pane_index} #{pane_title} "

tmux select-pane -t "$SESSION:work.0" -T "SIM"
tmux select-pane -t "$SESSION:work.2" -T "PARAVIEW"
tmux select-pane -t "$SESSION:work.1" -T "TEMPTRACK"

tmux select-pane -t "$SESSION:work.0"
tmux attach -t "$SESSION"
