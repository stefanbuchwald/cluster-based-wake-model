// Structured (transfinite) multi-block mesh for: cylinder in a channel
// Hex mesh via extrusion in z
//
// What you get:
// - Fully structured quad mesh in 2D (multi-block), then extruded to hex in 3D
// - Patch names: inlet, outlet, walls, cylinder, front, back
//
// Notes:
// - This is intentionally simple/robust: an inner square (size 2S) around the cylinder,
//   plus 12 outer blocks to reach the far channel boundaries.
// - Physical groups are created mostly via bounding boxes.
//
// (I cannot run gmsh here, so if gmsh complains about progression direction, simply flip the
//  endpoints of the corresponding Line(...) or set Pr=1.0.)

SetFactory("Built-in");

// ------------------------------------------------------------
// Zielgröße des Netzes
// ------------------------------------------------------------
// Setze hier direkt die gewünschte Gesamtzellzahl in Millionen:
//   1.0  -> ~0.98 Mio.
//   2.0  -> ~2.00 Mio.
//   3.0  -> ~3.00 Mio.
//
// Die Skalierung wirkt auf alle Transfinite-Divisionen in x/y,
// sodass du sauber zwischen grob / mittel / fein hochskalieren kannst.
TargetCellsM = 0.25;
Scale        = Sqrt(TargetCellsM);

// Basiskonfiguration (~0.98 Mio. Hexaeder bei TargetCellsM = 1.0)
Ntheta0 = 36;   // Punkte pro 45°-Bogen bzw. innerem Square-Segment
Nr0     = 24;   // Punkte radial Zylinder -> inneres Quadrat
NxUp0   = 15;   // Punkte in x von inlet bis x=-S
NxDn0   = 220;  // Punkte in x von x=+S bis outlet (mit Wake-Bias)
NyOut0  = 14;   // Punkte von y=±S bis zu den Kanalwänden
Nz0     = 52;   // Schichten in z
// -------------------------
// Parameters
// -------------------------
D      = 1.0;      // cylinder diameter
R      = D/2.0;    // cylinder radius

Lup    = 10.0*D;   // inlet length upstream of cylinder center
Ldown  = 25.0*D;   // outlet length downstream of cylinder center
H      = 10.0*D;   // full channel height

Lz     = 4.0*D;    // spanwise length
Nz     = Nz0;       // layers in z

S      = 2.5*D;    // half-size of inner square block (near-cylinder extent)

// -------------------------
// Mesh resolution
// -------------------------
// Rundung explizit, damit Gmsh ganzzahlige Transfinite-Werte bekommt
Ntheta = Floor(Ntheta0 * Scale + 0.5);
Nr     = Floor(Nr0     * Scale + 0.5);
NxUp   = Floor(NxUp0   * Scale + 0.5);
NxDn   = Floor(NxDn0   * Scale + 0.5);
NyOut  = Floor(NyOut0  * Scale + 0.5);

// Bias / Progression
// Wichtig:
// - Das frühere PrWake=1.03 bei NxDn=220 erzeugt am Nachlaufanfang extrem kleine
//   erste Zellen (~1e-3 D bei TargetCellsM=1) und war die Hauptursache für die
//   schlechten checkMesh-Metriken.
// - Daher werden die Progressionen an die globale Skalierung gekoppelt: mehr
//   Zellen => etwas schwächere Streckung, damit die kleinsten Zellen nicht
//   überproportional kollabieren.
PrRadBase  = 0.13;
PrUpBase   = 0.06;
PrWakeBase = 0.008;
PrWallBase = 0.12;

PrRad  = 1 + PrRadBase  / Scale; // radial: fein am Zylinder, gröber zum inner square
PrUp   = 1 + PrUpBase   / Scale; // upstream: fein Richtung Körper, gröber am inlet
PrWake = 1 + PrWakeBase / Scale; // downstream: fein im Nahnachlauf, moderat gröber zum outlet
PrWall = 1 + PrWallBase / Scale; // quer zur Strömung: Richtung Kanalrand gröber

// -------------------------
// Helper constants
// -------------------------
rt2 = Sqrt(2);

// -------------------------
// Points
// -------------------------
// Outer rectangle corner points
Point(1) = {-Lup,   -H/2, 0, 1.0};
Point(2) = { Ldown, -H/2, 0, 1.0};
Point(3) = { Ldown,  H/2, 0, 1.0};
Point(4) = {-Lup,    H/2, 0, 1.0};

// Outer boundary split points
// Inlet x=-Lup at y={-S,0,+S}
Point( 11) = {-Lup, -S,   0, 1.0};
Point( 12) = {-Lup,  0,   0, 1.0};
Point( 13) = {-Lup,  S,   0, 1.0};

// Outlet x=+Ldown at y={-S,0,+S}
Point( 21) = { Ldown, -S, 0, 1.0};
Point( 22) = { Ldown,  0, 0, 1.0};
Point( 23) = { Ldown,  S, 0, 1.0};

// Bottom y=-H/2 at x={-S,0,+S}
Point( 31) = {-S,  -H/2, 0, 1.0};
Point( 32) = { 0,  -H/2, 0, 1.0};
Point( 33) = { S,  -H/2, 0, 1.0};

// Top y=+H/2 at x={-S,0,+S}
Point( 41) = {-S,   H/2, 0, 1.0};
Point( 42) = { 0,   H/2, 0, 1.0};
Point( 43) = { S,   H/2, 0, 1.0};

// Inner square points (corners + midpoints on each side)
// Order around: (S,0)->(S,S)->(0,S)->(-S,S)->(-S,0)->(-S,-S)->(0,-S)->(S,-S)
Point(101) = { S,  0, 0, 1.0};
Point(102) = { S,  S, 0, 1.0};
Point(103) = { 0,  S, 0, 1.0};
Point(104) = {-S,  S, 0, 1.0};
Point(105) = {-S,  0, 0, 1.0};
Point(106) = {-S, -S, 0, 1.0};
Point(107) = { 0, -S, 0, 1.0};
Point(108) = { S, -S, 0, 1.0};

// Cylinder center + 8 points (every 45°)
Point(200) = { 0,  0, 0, 1.0};
Point(201) = { R,      0, 0, 1.0};
Point(202) = { R/rt2,  R/rt2, 0, 1.0};
Point(203) = { 0,      R, 0, 1.0};
Point(204) = {-R/rt2,  R/rt2, 0, 1.0};
Point(205) = {-R,      0, 0, 1.0};
Point(206) = {-R/rt2, -R/rt2, 0, 1.0};
Point(207) = { 0,     -R, 0, 1.0};
Point(208) = { R/rt2, -R/rt2, 0, 1.0};

// -------------------------
// Curves
// -------------------------
// Cylinder arcs (8x 45°)
Circle(301) = {201,200,202};
Circle(302) = {202,200,203};
Circle(303) = {203,200,204};
Circle(304) = {204,200,205};
Circle(305) = {205,200,206};
Circle(306) = {206,200,207};
Circle(307) = {207,200,208};
Circle(308) = {208,200,201};

// Inner-square perimeter (8 segments)
Line(401) = {101,102};
Line(402) = {102,103};
Line(403) = {103,104};
Line(404) = {104,105};
Line(405) = {105,106};
Line(406) = {106,107};
Line(407) = {107,108};
Line(408) = {108,101};

// Connections circle -> inner square (8 lines)
Line(501) = {201,101};
Line(502) = {202,102};
Line(503) = {203,103};
Line(504) = {204,104};
Line(505) = {205,105};
Line(506) = {206,106};
Line(507) = {207,107};
Line(508) = {208,108};

// Outer boundary (split)
// Bottom y=-H/2: 1 -> 31 -> 32 -> 33 -> 2
Line(601) = {1,31};
Line(602) = {31,32};
Line(603) = {32,33};
Line(604) = {33,2};

// Right x=+Ldown: 2 -> 21 -> 22 -> 23 -> 3
Line(611) = {2,21};
Line(612) = {21,22};
Line(613) = {22,23};
Line(614) = {23,3};

// Top y=+H/2: 3 -> 43 -> 42 -> 41 -> 4 (right->left)
Line(621) = {3,43};
Line(622) = {43,42};
Line(623) = {42,41};
Line(624) = {41,4};

// Left x=-Lup: 4 -> 13 -> 12 -> 11 -> 1 (top->bottom)
Line(631) = {4,13};
Line(632) = {13,12};
Line(633) = {12,11};
Line(634) = {11,1};

// Partition lines connecting inner square to outer boundary
// Horizontal y=+S: 13 -> 104   and   102 -> 23
Line(701) = {13,104};
Line(702) = {102,23};

// Horizontal y=-S: 11 -> 106 and 108 -> 21
Line(703) = {11,106};
Line(704) = {108,21};

// Vertical x=-S: 41 -> 104 and 106 -> 31
Line(705) = {41,104};
Line(706) = {106,31};

// Vertical x=+S: 43 -> 102 and 108 -> 33
Line(707) = {43,102};
Line(708) = {108,33};

// Vertical x=0: 42 -> 103 and 107 -> 32
Line(709) = {42,103};
Line(710) = {107,32};

// Horizontal y=0: 12 -> 105 and 101 -> 22
Line(711) = {12,105};
Line(712) = {101,22};

// -------------------------
// Surfaces (8 inner + 12 outer)
// -------------------------

// --- Inner (between cylinder and inner square), 8 quads ---
Line Loop(801) = { 301,  502, -401, -501};
Plane Surface(901) = {801};
Line Loop(802) = { 302,  503, -402, -502};
Plane Surface(902) = {802};
Line Loop(803) = { 303,  504, -403, -503};
Plane Surface(903) = {803};
Line Loop(804) = { 304,  505, -404, -504};
Plane Surface(904) = {804};
Line Loop(805) = { 305,  506, -405, -505};
Plane Surface(905) = {805};
Line Loop(806) = { 306,  507, -406, -506};
Plane Surface(906) = {806};
Line Loop(807) = { 307,  508, -407, -507};
Plane Surface(907) = {807};
Line Loop(808) = { 308,  501, -408, -508};
Plane Surface(908) = {808};

// --- Outer 12 blocks ---
// Top-left corner: vertices 4 -> 41 -> 104 -> 13
Line Loop(811) = {-624, 705, -701, -631};
Plane Surface(911) = {811};

// Top-mid-left: vertices 41 -> 42 -> 103 -> 104
Line Loop(812) = {-623, 709, 403, -705};
Plane Surface(912) = {812};

// Top-mid-right: vertices 42 -> 43 -> 102 -> 103
Line Loop(813) = {-622, 707, 402, -709};
Plane Surface(913) = {813};

// Top-right corner: vertices 43 -> 3 -> 23 -> 102
Line Loop(814) = {-621, -614, -702, -707};
Plane Surface(914) = {814};

// Bottom-left corner: vertices 1 -> 31 -> 106 -> 11
Line Loop(815) = {601, -706, -703, 634};
Plane Surface(915) = {815};

// Bottom-mid-left: vertices 31 -> 32 -> 107 -> 106
Line Loop(816) = {602, -710, -406, 706};
Plane Surface(916) = {816};

// Bottom-mid-right: vertices 32 -> 33 -> 108 -> 107
Line Loop(817) = {603, -708, -407, 710};
Plane Surface(917) = {817};

// Bottom-right corner: vertices 33 -> 2 -> 21 -> 108
Line Loop(818) = {604, 611, -704, 708};
Plane Surface(918) = {818};

// Left-middle-top: vertices 12 -> 13 -> 104 -> 105
Line Loop(819) = {-632, 701, 404, -711};
Plane Surface(919) = {819};

// Left-middle-bottom: vertices 11 -> 12 -> 105 -> 106
Line Loop(820) = {-633, 711, 405, -703};
Plane Surface(920) = {820};

// Right-middle-top: vertices 101 -> 102 -> 23 -> 22
Line Loop(821) = {401, 702, -613, -712};
Plane Surface(921) = {821};

// Right-middle-bottom: vertices 108 -> 101 -> 22 -> 21
Line Loop(822) = {408, 712, -612, -704};
Plane Surface(922) = {822};

// -------------------------
// Transfinite definitions
// -------------------------

// Cylinder arcs
Transfinite Curve {301,302,303,304,305,306,307,308} = Ntheta;

// Inner square perimeter (8 segments) -> same as arcs
Transfinite Curve {401,402,403,404,405,406,407,408} = Ntheta;

// Connections cylinder->square
Transfinite Curve {501,502,503,504,505,506,507,508} = Nr Using Progression PrRad;

// Outer boundary segments
// x-Richtung: außen grober, in Körpernähe / Nachlaufanfang feiner
Transfinite Curve {601} = NxUp Using Progression 1/PrUp;
Transfinite Curve {602} = Ntheta;
Transfinite Curve {603} = Ntheta;
Transfinite Curve {604} = NxDn Using Progression PrWake;

Transfinite Curve {621} = NxDn Using Progression 1/PrWake;
Transfinite Curve {622} = Ntheta;
Transfinite Curve {623} = Ntheta;
Transfinite Curve {624} = NxUp Using Progression PrUp;

// y-Richtung: an den Kanalrändern grober, Richtung Kernzone feiner
Transfinite Curve {631} = NyOut Using Progression 1/PrWall;
Transfinite Curve {632} = Ntheta;
Transfinite Curve {633} = Ntheta;
Transfinite Curve {634} = NyOut Using Progression PrWall;

Transfinite Curve {611} = NyOut Using Progression 1/PrWall;
Transfinite Curve {612} = Ntheta;
Transfinite Curve {613} = Ntheta;
Transfinite Curve {614} = NyOut Using Progression PrWall;

// Partitions to connect inner square to outer boundary
Transfinite Curve {701,703,711} = NxUp Using Progression 1/PrUp;
Transfinite Curve {702,704,712} = NxDn Using Progression PrWake;

Transfinite Curve {705,707,709} = NyOut Using Progression 1/PrWall;
Transfinite Curve {706,708,710} = NyOut Using Progression PrWall;

// Apply transfinite surfaces (corner points specified for robustness)
Transfinite Surface {901} = {201,202,102,101};
Transfinite Surface {902} = {202,203,103,102};
Transfinite Surface {903} = {203,204,104,103};
Transfinite Surface {904} = {204,205,105,104};
Transfinite Surface {905} = {205,206,106,105};
Transfinite Surface {906} = {206,207,107,106};
Transfinite Surface {907} = {207,208,108,107};
Transfinite Surface {908} = {208,201,101,108};

Transfinite Surface {911} = {4,41,104,13};
Transfinite Surface {912} = {41,42,103,104};
Transfinite Surface {913} = {42,43,102,103};
Transfinite Surface {914} = {43,3,23,102};

Transfinite Surface {915} = {1,31,106,11};
Transfinite Surface {916} = {31,32,107,106};
Transfinite Surface {917} = {32,33,108,107};
Transfinite Surface {918} = {33,2,21,108};

Transfinite Surface {919} = {12,13,104,105};
Transfinite Surface {920} = {11,12,105,106};

Transfinite Surface {921} = {101,102,23,22};
Transfinite Surface {922} = {108,101,22,21};

// Recombine all 2D surfaces to quads
Recombine Surface {901:908, 911:922};

// -------------------------
// Extrude to 3D (hex)
// -------------------------
// Note: In Gmsh >= 4.11 the old 'Compound Surface(id) = {...};' syntax is deprecated.
// We extrude all transfinite quad blocks; they share curves so the 3D mesh is conformal.
//Extrude {0,0,Lz} {
//  Surface{901:908, 911:922};
//  Layers{Nz};
//  Recombine;
//};

// -------------------------
// Extrude to 3D (hex) - robust for Gmsh 4.12
// -------------------------

baseSurfs[] = {901,902,903,904,905,906,907,908, 911,912,913,914,915,916,917,918,919,920,921,922};

// (optional) ensure all are quads before extrusion
Recombine Surface{baseSurfs[]};

// IMPORTANT: put everything in ONE statement (no range syntax)
out[] = Extrude {0, 0, Lz} {
  Surface{baseSurfs[]};
  Layers{Nz};
  Recombine;
};

// -------------------------
// Physical groups (OpenFOAM)
// -------------------------
xtol = 1e-6;
ztol = 1e-6;
rtol = 1e-3*D;

// Volume: entire domain (select by bounding box)
Physical Volume("fluid") =
{ Volume In BoundingBox{-Lup-1, -H/2-1, -ztol, Ldown+1, H/2+1, Lz+ztol} };

// z-normal planes at z=0 and z=Lz (mini-LES spanwise coupling happens in OpenFOAM via cyclic)
Physical Surface("front") =
{ Surface In BoundingBox{-Lup-1, -H/2-1, -ztol, Ldown+1, H/2+1,  ztol} };
Physical Surface("back")  =
{ Surface In BoundingBox{-Lup-1, -H/2-1, Lz-ztol, Ldown+1, H/2+1, Lz+ztol} };

// inlet/outlet (x = -Lup / +Ldown)
Physical Surface("inlet")  =
{ Surface In BoundingBox{-Lup-xtol, -H/2-xtol, -ztol, -Lup+xtol,  H/2+xtol, Lz+ztol} };
Physical Surface("outlet") =
{ Surface In BoundingBox{ Ldown-xtol, -H/2-xtol, -ztol,  Ldown+xtol, H/2+xtol, Lz+ztol} };

// channel walls (y = ±H/2)
Physical Surface("walls") =
{
  Surface In BoundingBox{-Lup-1,  H/2-xtol, -ztol, Ldown+1, H/2+xtol, Lz+ztol},
  Surface In BoundingBox{-Lup-1, -H/2-xtol, -ztol, Ldown+1,-H/2+xtol, Lz+ztol}
};

// cylinder wall (tight box around radius R)
Physical Surface("cylinder") =
{ Surface In BoundingBox{-R-rtol, -R-rtol, -ztol, R+rtol, R+rtol, Lz+ztol} };

// -------------------------
// Geschätzte Zellzahl
// -------------------------
// Transfinite Curve zählt Punkte, daher Zellzahl = (N-1) je Richtung.
QuadCells =
  8 * (Ntheta - 1) * (Nr    - 1)
+ 2 * (NxUp  - 1) * (NyOut - 1)
+ 2 * (NxDn  - 1) * (NyOut - 1)
+ 4 * (Ntheta- 1) * (NyOut - 1)
+ 2 * (Ntheta- 1) * (NxUp  - 1)
+ 2 * (Ntheta- 1) * (NxDn  - 1);

HexCells = QuadCells * Nz;

Printf("TargetCellsM = %g", TargetCellsM);
Printf("Approx. hex cells = %g", HexCells);
Printf("PrRad = %g", PrRad);
Printf("PrUp = %g", PrUp);
Printf("PrWake = %g", PrWake);
Printf("PrWall = %g", PrWall);

// -------------------------
// Output options
// -------------------------
Mesh.MshFileVersion = 2.2;
