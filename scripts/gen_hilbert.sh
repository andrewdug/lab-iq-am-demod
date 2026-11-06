#!/bin/bash
# Generates Qucs-S schematics of FIR Hilbert filters (MicroCap ex5 analogue):
#   H  - order 10 (N = 11, n0 = 6), rectangular window      -> node y,   reference xd   (5 Dt)
#   W  - order 10, Hamming window                            -> node yw,  reference xd   (5 Dt)
#   O  - order 14 (N = 15, n0 = 8), Hamming window (p. 4, final filter) -> node y14, reference xd14 (7 Dt)
#   G  - order 18 (N = 19, n0 = 10), Hamming window (margin check)      -> node y18, reference xd18 (9 Dt)
# Usage: gen_hilbert.sh <mode: AC|TR> <out.sch>
set -euo pipefail
mode=${1:?mode AC|TR}; out=${2:?output .sch}; name=$(basename "$out" .sch)
case $mode in AC|TR) ;; *) echo "mode must be AC or TR" >&2; exit 1;; esac

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
comps=$tmp/comps; wires=$tmp/wires; diags=$tmp/diags
: > "$comps"; : > "$wires"; : > "$diags"

# Layout in the style of MicroCap ex5 (all connections are real wires; Qucs joins wires only at end points):
#   delay line D1..D9 (Ideal/VDelay, 2Dt each) on y = 70; tap bus of node x(2j) runs down at x = 350 + 200 j;
#   filter rows: G at y = 500, O at 750, W at 1000, H at 1250 (longest first, so extra taps do not cross the other rows);
#   tap j: gain block (VCVS, T = 0) at x = 400 + 200 j, adder (Ideal/VSum) A_j at (500 + 200 j, row + 70);
#   reference delays D0 (5 Dt), D14 (7 Dt) and D18 (9 Dt) at the bottom, fed from the x0 bus.
# VCVS pins: in+ (x-30, y-30), in- (x-30, y+30), out+ (x+30, y-30), out- (x+30, y+30).
# VDelay pins: in (x-30, y), out (x+30, y).  VSum pins: in1 (x, y-30), in2 (x, y+30), out (x+30, y).
ROWS=(500 750 1000 1250)
ROW_LAST_TAP=(9 7 5 5)   # last tap index j used by each row

# wire X1 Y1 X2 Y2 [LABEL]
wire() {
  local l=${5:-}
  if [ -n "$l" ]; then
    echo "  <$1 $2 $3 $4 \"$l\" $(( ($1 + $3) / 2 - 10 )) $(( $2 - 30 )) 0 \"\">" >> "$wires"
  else
    echo "  <$1 $2 $3 $4 \"\" 0 0 0 \"\">" >> "$wires"
  fi
}
gnd() { echo "  <GND * 1 $1 $2 0 0 0 0>" >> "$comps"; }
# text X Y SIZE "TEXT" - annotation (Paintings)
text() { echo "  <Text $1 $2 $3 #00007f 0 \"$4\">" >> "$diags"; }
# gain NAME X Y K: ideal amplifier (VCVS without delay), in- and out- grounded
gain() {
  echo "  <VCVS $1 1 $2 $3 -26 34 0 0 \"$4\" 1 \"0\" 0>" >> "$comps"
  gnd $(($2 - 30)) $(($3 + 30)); gnd $(($2 + 30)) $(($3 + 30))
}
# delay NAME X Y TDEL: Ideal/VDelay block
delay() { echo "  <Lib $1 1 $2 $3 -20 24 0 0 \"Ideal\" 0 \"VDelay\" 0 \"$4\" 1>" >> "$comps"; }
# adder NAME X Y: Ideal/VSum block
adder() { echo "  <Lib $1 1 $2 $3 22 -36 0 0 \"Ideal\" 0 \"VSum\" 0>" >> "$comps"; }
# rload NAME X Y: vertical 1 MOhm, top pin (X, Y-30), bottom grounded
rload() { echo "  <R $1 1 $2 $3 15 -26 0 1 \"1 MOhm\" 1 \"26.85\" 0 \"0.0\" 0 \"0.0\" 0 \"26.85\" 0 \"european\" 0>" >> "$comps"
  gnd "$2" $(($3 + 30)); }
# chain PREFIX ROW ORDER COEF OUTNODE LOADNAME: taps x0, x2, ..., x<ORDER>;
# y = K0 x0 + K1 x2 + ... accumulated by the adder chain A1 ... A<nt-1>
chain() {
  local p=$1 row=$2 order=$3 cf=$4 outn=$5 ln=$6 nt=$(( $3 / 2 + 1 )) j x ax
  local ay=$((row + 70)) yb=$((row + 130))
  for ((j = 0; j < nt; j++)); do
    x=$((400 + 200 * j)); ax=$((x + 100))
    gain "$p$j" $x "$row" "$cf$j"
    wire $((x - 50)) $((row - 30)) $((x - 30)) $((row - 30))          # tap bus -> gain input
    wire $((x + 30)) $((row - 30)) $ax $((row - 30))                   # gain output -> right
    if [ $j -eq 0 ]; then                                              # K0 x0 -> lower input of A1
      wire $ax $((row - 30)) $ax $yb
      wire $ax $yb $((ax + 200)) $yb
      wire $((ax + 200)) $yb $((ax + 200)) $((ay + 30))
    else
      adder "A$p$j" $ax $ay
      wire $ax $((row - 30)) $ax $((ay - 30))                          # K_j x(2j) -> upper input of A_j
      if [ $j -lt $((nt - 1)) ]; then                                  # A_j output -> lower input of A_(j+1)
        wire $((ax + 30)) $ay $((ax + 30)) $yb
        wire $((ax + 30)) $yb $((ax + 200)) $yb
        wire $((ax + 200)) $yb $((ax + 200)) $((ay + 30))
      fi
    fi
  done
  ax=$((400 + 200 * (nt - 1) + 100))
  wire $((ax + 30)) $ay $((ax + 100)) $ay "$outn"
  rload "$ln" $((ax + 100)) $((ay + 30))
}

# --- source ---
if [ "$mode" = TR ]; then
  echo '  <AM_Mod V1 1 100 200 18 -26 0 0 "U0" 1 "f0" 1 "0" 1 "m" 1>' >> "$comps"
  echo '  <Vac Vm 1 -20 230 18 -26 0 1 "1" 1 "Fm" 1 "0" 0 "0" 0 "0" 0 "0" 0>' >> "$comps"
  echo '  <-20 200 70 200 "" 0 0 0 "">' >> "$wires"
  echo '  <-20 260 100 260 "" 0 0 0 "">' >> "$wires"
  echo '  <100 230 100 260 "" 0 0 0 "">' >> "$wires"
  gnd 100 260
else
  echo '  <Vac V1 1 100 200 18 -26 0 1 "1 V" 1 "1 GHz" 0 "0" 0 "0" 0 "0" 0 "0" 0>' >> "$comps"
  gnd 100 230
fi
wire 100 170 100 70

# --- delay line D1..D9: nine 2*Dt blocks, nodes x0, x2, ..., x18 = s(t - k*Dt) ---
wire 100 70 350 70 x0
for i in 1 2 3 4 5 6 7 8 9; do
  x=$((200 + 200 * i))
  delay D$i $x 70 T2
  wire $((x - 50)) 70 $((x - 30)) 70                                   # node x(2i-2) -> in
  wire $((x + 30)) 70 $((x + 150)) 70 "x$((2 * i))"                     # out -> node x(2i)
done

# --- tap buses: node x(2j) down to the filter rows (split at every row so the joints are wire ends) ---
for j in 0 1 2 3 4 5 6 7 8 9; do
  bx=$((350 + 200 * j)); y0=70
  for r in "${!ROWS[@]}"; do
    row=${ROWS[$r]}
    [ $j -gt "${ROW_LAST_TAP[$r]}" ] && continue
    wire $bx $y0 $bx $((row - 30)); y0=$((row - 30))
  done
done

# --- filters ---
chain G 500 18 g y18 Ry18
chain O 750 14 o y14 Ry14
chain W 1000 10 cw yw Ryw
chain H 1250 10 c y Ry

# --- delay line of the source signal s'(t) (p. 5): D0 = 5 Dt -> xd, D0 + D14 = 5 Dt + 2 Dt -> xd14, D18 = 9 Dt -> xd18 ---
# A single VDelay with Tdel = 7 Dt gives isolated wrong points in qucsator TR (up to 0.25 V at 119 of 41281 points),
# while 2 Dt, 5 Dt and 9 Dt are clean to 1e-12 V, so 7 Dt is built as a cascade of the two clean delays.
wire 350 1220 350 1550; wire 350 1550 350 1790
delay D0 450 1550 T5;  wire 350 1550 420 1550; wire 480 1550 600 1550 xd;   rload Rxd 600 1580
delay D14 700 1550 T2; wire 600 1550 670 1550; wire 730 1550 850 1550 xd14; rload Rxd14 850 1580
delay D18 450 1790 T9; wire 350 1790 420 1790; wire 480 1790 600 1790 xd18; rload Rxd18 600 1820

# --- annotations ---
text -300 -110 14 "ЛИНИЯ ЗАДЕРЖКИ ФИЛЬТРА: D1...D9 - блоки задержки VDelay (библиотека Ideal), Tdel = T2 = 2Dt = 1,163 нс.\nУзлы x0, x2, ..., x18 - отсчёты s(t - k Dt); вертикальные шины подают их на отводы фильтров."
text -300 430 14 "ПОРЯДОК 18, окно Хэмминга (проверка запаса):\nусилители G0...G9, K = g0...g9,\nсумматоры AG1...AG9 -> y18"
text -300 680 14 "ПОРЯДОК 14, окно Хэмминга (п. 4, итоговый):\nусилители O0...O7, K = o0...o7,\nсумматоры AO1...AO7 -> y14"
text -300 930 14 "ПОРЯДОК 10, окно Хэмминга:\nусилители W0...W5, K = cw0...cw5,\nсумматоры AW1...AW5 -> yw"
text -300 1180 14 "ПОРЯДОК 10, без окна:\nусилители H0...H5, K = c0...c5,\nсумматоры AH1...AH5 -> y"
text -300 1510 14 "ЛИНИЯ ЗАДЕРЖКИ ИСХОДНОГО\nСИГНАЛА s'(t) (п. 5):\nD0: 5Dt = 2,907 нс -> xd (к порядку 10)\nD0 + D14: 5Dt + 2Dt = 4,070 нс -> xd14 (к порядку 14)\nD18: 9Dt = 5,233 нс -> xd18 (к порядку 18)"
text 1000 1530 14 "Усилитель с коэффициентом K - ИНУН (VCVS) без задержки; VDelay и VSum внутри тоже ИНУН.\nR = 1 МОм - нагрузки, только для вывода напряжений y, yw, y14, y18, xd, xd14, xd18.\nДискретизация с шагом Dt и п. 5-8 (z(n), годограф, спектры, I/Q, демодуляция):\nPython-скрипт scripts/hilbert_iq.py, запускать после моделирования AC и TR."
# order-18 Hamming coefficients: g_j = 2/(pi k) * (0.54 - 0.46 cos(2 pi n / 18)), n = 2j, k = n - 9
geq=""
for j in 0 1 2 3 4 5 6 7 8 9; do
  geq="$geq \"g$j=2/(pi*($((2 * j - 9))))*(0.54-0.46*cos(2*pi*$((2 * j))/18))\" 1"
done
# order-14 Hamming coefficients: o_j = 2/(pi k) * (0.54 - 0.46 cos(2 pi n / 14)), n = 2j, k = n - 7
oeq=""
for j in 0 1 2 3 4 5 6 7; do
  oeq="$oeq \"o$j=2/(pi*($((2 * j - 7))))*(0.54-0.46*cos(2*pi*$((2 * j))/14))\" 1"
done
onum=$(awk 'BEGIN{pi=atan2(0,-1); for(j=0;j<8;j++){n=2*j; k=n-7; printf "%s%.4f", (j?" ":""), 2/(pi*k)*(0.54-0.46*cos(2*pi*n/14))}}' | tr . ,)
gnum=$(awk 'BEGIN{pi=atan2(0,-1); for(j=0;j<10;j++){n=2*j; k=n-9; printf "%s%.4f", (j?" ":""), 2/(pi*k)*(0.54-0.46*cos(2*pi*n/18))}}' | tr . ,)

cat >> "$comps" <<EOF
  <Eqn Eqn1 1 -300 -750 -37 18 0 0 "N=1" 1 "f0=(400+30*N)*1e6" 1 "Fm=3*N*1e6" 1 "U0=1" 1 "m=1" 1 "fs=4*f0" 1 "Dt=1/fs" 1 "T2=2*Dt" 1 "T5=5*Dt" 1 "T7=7*Dt" 1 "T9=9*Dt" 1 "yes" 0>
  <Eqn EqnH 1 200 -750 -37 18 0 0 "c0=-2/(5*pi)" 1 "c1=-2/(3*pi)" 1 "c2=-2/pi" 1 "c3=2/pi" 1 "c4=2/(3*pi)" 1 "c5=2/(5*pi)" 1 "w5=0.54-0.46*cos(2*pi*0/10)" 1 "w3=0.54-0.46*cos(2*pi*2/10)" 1 "w1=0.54-0.46*cos(2*pi*4/10)" 1 "cw0=c0*w5" 1 "cw1=c1*w3" 1 "cw2=c2*w1" 1 "cw3=c3*w1" 1 "cw4=c4*w3" 1 "cw5=c5*w5" 1 "yes" 0>
  <Eqn EqnG 1 700 -750 -37 18 0 0$geq "yes" 0>
  <Eqn EqnO 1 1250 -750 -37 18 0 0$oeq "yes" 0>
EOF

if [ "$mode" = TR ]; then
  # 41281 points: step Dt/12, so the samples t = n*Dt lie exactly on the output grid
  cat >> "$comps" <<'EOF'
  <.TR TR1 1 -300 1960 0 78 0 0 "lin" 1 "0" 1 "2e-6" 1 "41281" 1 "Trapezoidal" 0 "2" 0 "1 ps" 0 "1e-16" 0 "150" 0 "0.001" 0 "1 pA" 0 "1 uV" 0 "26.85" 0 "1e-3" 0 "1e-6" 0 "1" 0 "CroutLU" 0 "no" 0 "yes" 0 "10 ps" 0>
  <Eqn EqnE 1 1800 -750 -37 18 0 0 "Uenv=U0*(1+m*sin(2*pi*Fm*(time-T5)))" 1 "Uenv14=U0*(1+m*sin(2*pi*Fm*(time-T7)))" 1 "Uenv18=U0*(1+m*sin(2*pi*Fm*(time-T9)))" 1 "Ur=sqrt(xd.Vt^2+y.Vt^2)" 1 "Uw=sqrt(xd.Vt^2+yw.Vt^2)" 1 "Uo=sqrt(xd14.Vt^2+y14.Vt^2)" 1 "Ug=sqrt(xd18.Vt^2+y18.Vt^2)" 1 "yes" 0>
EOF
  cat >> "$diags" <<'EOF'
  <Rect 3000 -300 1000 300 3 #c0c0c0 1 00 1 0 1e-07 1e-06 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "t, с" "В" "">
	<"x0.Vt" #0000ff 0 3 0 0 0>
  </Rect>
  <Rect 3000 200 1000 300 3 #c0c0c0 1 00 0 1e-06 2e-09 1.01e-06 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "t, с" "В" "">
	<"x0.Vt" #008000 0 3 0 0 0>
	<"xd.Vt" #000000 1 3 0 0 0>
	<"y.Vt" #ff0000 1 3 0 0 0>
	<"yw.Vt" #0000ff 1 3 0 1 0>
  </Rect>
  <Rect 3000 700 1000 300 3 #c0c0c0 1 00 0 1e-06 2e-09 1.01e-06 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "t, с" "В" "">
	<"x0.Vt" #008000 0 3 0 0 0>
	<"xd18.Vt" #000000 1 3 0 0 0>
	<"y18.Vt" #0000ff 1 3 0 0 0>
  </Rect>
  <Rect 3000 1200 1000 300 3 #c0c0c0 1 00 1 0 1e-07 1e-06 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "t, с" "В" "">
	<"Uenv" #000000 2 3 0 0 0>
	<"Ur" #ff0000 0 3 0 0 0>
	<"Uw" #0000ff 0 3 0 0 0>
  </Rect>
  <Rect 3000 1700 1000 300 3 #c0c0c0 1 00 1 0 1e-07 1e-06 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "t, с" "В" "">
	<"Uenv18" #000000 2 3 0 0 0>
	<"Ug" #0000ff 0 3 0 0 0>
  </Rect>
  <Rect 3000 2200 1000 300 3 #c0c0c0 1 00 0 1e-06 2e-09 1.01e-06 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "t, с" "В" "">
	<"x0.Vt" #008000 0 3 0 0 0>
	<"xd14.Vt" #000000 1 3 0 0 0>
	<"y14.Vt" #ff6600 1 3 0 0 0>
  </Rect>
  <Rect 3000 2700 1000 300 3 #c0c0c0 1 00 1 0 1e-07 1e-06 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "t, с" "В" "">
	<"Uenv14" #000000 2 3 0 0 0>
	<"Uo" #ff6600 0 3 0 0 0>
  </Rect>
  <Text 3000 1860 14 #00007f 0 "Порядок 14, окно Хэмминга (итоговый): s(t) (зелёный), s'(t) = s(t - 7Dt) = xd14 (чёрный), сопряжённый y14 (оранжевый)">
  <Text 3000 2360 14 #00007f 0 "Порядок 14: огибающая sqrt(xd14^2 + y14^2) (оранжевый) и модулирующее U0(1 + m sin 2pi Fм (t - 7Dt)) (чёрный)">
  <Text 3000 -640 14 #00007f 0 "s(t) = x0: входной АМ-сигнал, f0 = 430 МГц, Fм = 3 МГц, m = 1, U0 = 1 В">
  <Text 3000 -140 14 #00007f 0 "Порядок 10: s(t) (зелёный), s'(t) = s(t - 5Dt) = xd (чёрный), сопряжённый y - без окна (красный), yw - окно Хэмминга (синий пунктир)">
  <Text 3000 360 14 #00007f 0 "Порядок 18, окно Хэмминга: s(t) (зелёный), s'(t) = s(t - 9Dt) = xd18 (чёрный), сопряжённый y18 (синий)">
  <Text 3000 860 14 #00007f 0 "Порядок 10: огибающая sqrt(xd^2 + y^2) без окна (красный), sqrt(xd^2 + yw^2) с окном (синий), модулирующее U0(1 + m sin 2pi Fм (t - 5Dt)) (чёрный)">
  <Text 3000 1360 14 #00007f 0 "Порядок 18: огибающая sqrt(xd18^2 + y18^2) (синий) и модулирующее U0(1 + m sin 2pi Fм (t - 9Dt)) (чёрный)">
EOF
else
  cat >> "$comps" <<'EOF'
  <.AC AC1 1 -300 1960 0 78 0 0 "lin" 1 "1 kHz" 1 "860 MHz" 1 "8601" 1 "no" 0>
  <Eqn EqnA 1 1800 -750 -37 18 0 0 "Kh=mag(y.v)" 1 "Kw=mag(yw.v)" 1 "Ko=mag(y14.v)" 1 "Kg=mag(y18.v)" 1 "Kd=mag(xd.v)" 1 "Fh=unwrap(angle(y.v))" 1 "Fw=unwrap(angle(yw.v))" 1 "Fo=unwrap(angle(y14.v))" 1 "Fg=unwrap(angle(y18.v))" 1 "Fd=unwrap(angle(xd.v))" 1 "dFh=angle(y.v/xd.v)*180/pi" 1 "dFw=angle(yw.v/xd.v)*180/pi" 1 "dFo=angle(y14.v/xd14.v)*180/pi" 1 "dFg=angle(y18.v/xd18.v)*180/pi" 1 "tgw=-diff(Fw,acfrequency)/(2*pi)" 1 "tgo=-diff(Fo,acfrequency)/(2*pi)" 1 "tgg=-diff(Fg,acfrequency)/(2*pi)" 1 "yes" 0>
  <Eqn EqnF0 1 2350 -750 -37 18 0 0 "If0=4300" 1 "Kh_f0=Kh[If0]" 1 "Kw_f0=Kw[If0]" 1 "Ko_f0=Ko[If0]" 1 "Kg_f0=Kg[If0]" 1 "dKh_f0=(Kh[If0]-Kd[If0])*100" 1 "dKw_f0=(Kw[If0]-Kd[If0])*100" 1 "dKo_f0=(Ko[If0]-Kd[If0])*100" 1 "dKg_f0=(Kg[If0]-Kd[If0])*100" 1 "dFh_f0=dFh[If0]" 1 "dFw_f0=dFw[If0]" 1 "dFo_f0=dFo[If0]" 1 "dFg_f0=dFg[If0]" 1 "t0w_ns=tgw[If0]*1e9" 1 "t0o_ns=tgo[If0]*1e9" 1 "t0g_ns=tgg[If0]*1e9" 1 "T5_ns=T5*1e9" 1 "T7_ns=T7*1e9" 1 "T9_ns=T9*1e9" 1 "yes" 0>
EOF
  cat >> "$diags" <<'EOF'
  <Rect 3000 -300 1000 350 3 #c0c0c0 1 00 0 0 1e+08 8.6e+08 0 0 0.2 1.4 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "|K|" "">
	<"Kd" #000000 1 3 0 1 0>
	<"Kh" #ff0000 2 3 0 0 0>
	<"Kw" #0000ff 2 3 0 0 0>
	<"Ko" #ff6600 2 3 0 0 0>
	<"Kg" #008000 2 3 0 0 0>
  </Rect>
  <Rect 3000 250 1000 350 3 #c0c0c0 1 00 0 0 1e+08 8.6e+08 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "рад" "">
	<"Fd" #000000 1 3 0 1 0>
	<"Fh" #ff0000 2 3 0 0 0>
	<"Fw" #0000ff 1 3 0 1 0>
	<"Fo" #ff6600 2 3 0 0 0>
	<"Fg" #008000 1 3 0 0 0>
  </Rect>
  <Rect 3000 800 1000 350 3 #c0c0c0 1 00 0 0 1e+08 8.5e+08 0 -180 30 0 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "град" "">
	<"dFh" #ff0000 2 3 0 0 0>
	<"dFw" #0000ff 1 3 0 1 0>
	<"dFo" #ff6600 2 3 0 0 0>
	<"dFg" #008000 1 3 0 0 0>
  </Rect>
  <Rect 3000 1350 1000 350 3 #c0c0c0 1 00 0 1e+08 1e+08 7.6e+08 0 0 1e-09 6e-09 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "с" "">
	<"tgw" #0000ff 1 3 0 0 0>
	<"tgo" #ff6600 1 3 0 0 0>
	<"tgg" #008000 1 3 0 0 0>
  </Rect>
  <Tab 3000 1850 1500 120 3 #c0c0c0 1 00 1 0 1 1 1 0 1 1 1 0 1 1 315 0 225 0 0 0 "" "" "">
	<"Kh_f0" #0000ff 0 4 0 0 0>
	<"Kw_f0" #0000ff 0 4 0 0 0>
	<"Ko_f0" #0000ff 0 4 0 0 0>
	<"Kg_f0" #0000ff 0 4 0 0 0>
	<"dKh_f0" #0000ff 0 3 0 0 0>
	<"dKw_f0" #0000ff 0 3 0 0 0>
	<"dKo_f0" #0000ff 0 3 0 0 0>
	<"dKg_f0" #0000ff 0 3 0 0 0>
	<"dFh_f0" #0000ff 0 3 0 0 0>
	<"dFw_f0" #0000ff 0 3 0 0 0>
	<"dFo_f0" #0000ff 0 3 0 0 0>
	<"dFg_f0" #0000ff 0 3 0 0 0>
	<"t0w_ns" #0000ff 0 4 0 0 0>
	<"T5_ns" #0000ff 0 4 0 0 0>
	<"t0o_ns" #0000ff 0 4 0 0 0>
	<"T7_ns" #0000ff 0 4 0 0 0>
	<"t0g_ns" #0000ff 0 4 0 0 0>
	<"T9_ns" #0000ff 0 4 0 0 0>
  </Tab>
  <Text 3000 -690 14 #00007f 0 "АЧХ: задержка (чёрный пунктир), порядок 10 без окна (красный), порядок 10 с окном Хэмминга (синий), порядок 14 с окном (оранжевый), порядок 18 с окном (зелёный); f0 = fs/4 = 430 МГц">
  <Text 3000 -140 14 #00007f 0 "ФЧХ (развёрнутая), рад: задержка 5Dt (чёрный), порядок 10 без окна (красный), с окном (синий), порядок 14 с окном (оранжевый), порядок 18 с окном (зелёный)">
  <Text 3000 410 14 #00007f 0 "Разность фаз выхода фильтра и своего элемента задержки (5Dt для порядка 10, 7Dt для порядка 14, 9Dt для порядка 18), град">
  <Text 3000 960 14 #00007f 0 "Групповая задержка tгр = -dФ/dw: порядок 10 (синий, теория 5Dt = 2,907 нс), порядок 14 (оранжевый, 7Dt = 4,070 нс), порядок 18 (зелёный, 9Dt = 5,233 нс)">
  <Text 3000 1510 14 #00007f 0 "На частоте f0 = 430 МГц: |K|, отклонение |K| от 1 (%), разность фаз (град), задержка по ФЧХ и теоретическая (нс)">
EOF
fi
cat >> "$diags" <<EOF
  <Text -300 -980 14 #00007f 0 "Цифровые КИХ-фильтры Гильберта (аналог модели MicroCap ex5), N = 1:\nfs = 4 f0 = 1720 МГц (f0 - в середине области равномерности АЧХ), Dt = 1/fs = 0,5814 нс; линия задержки 9 x 2Dt.\nПорядок 10 (11 коэфф., n0 = 6): h[k] = 2/(pi k), k = -5..5 нечётные: -0,1273 -0,2122 -0,6366 0 0,6366 0,2122 0,1273; задержка t0 = 5Dt = 2,907 нс.\nОкно Хэмминга w[n] = 0,54 - 0,46 cos(2 pi n/(N-1)), N = 11: 0,08 0,398 0,912 1 0,912 0,398 0,08.\nПорядок 14 (15 коэфф., n0 = 8) с окном Хэмминга, k = -7..7 нечётные: $onum; t0 = 7Dt = 4,070 нс.\nПорядок 18 (19 коэфф., n0 = 10) с окном Хэмминга, k = -9..9 нечётные: $gnum; t0 = 9Dt = 5,233 нс.\nЗадержки - блоки VDelay, умножение на коэффициент - ИНУН, сложение - цепочка VSum.">
EOF

{
cat <<EOF
<Qucs Schematic 25.2.99>
<Properties>
  <View=-400,-1000,4100,3200,0.5,0,0>
  <Grid=10,10,1>
  <DataSet=$name.dat>
  <DataDisplay=$name.dpl>
  <OpenDisplay=0>
  <Script=$name.m>
  <RunScript=0>
  <showFrame=0>
  <FrameText0=Название>
  <FrameText1=Чертил:>
  <FrameText2=Дата:>
  <FrameText3=Версия:>
</Properties>
<Symbol>
</Symbol>
<Components>
EOF
cat "$comps"
echo "</Components>"; echo "<Wires>"; cat "$wires"; echo "</Wires>"
echo "<Diagrams>"; grep -v '^  <Text' "$diags"; echo "</Diagrams>"; echo "<Paintings>"; grep '^  <Text' "$diags"; echo "</Paintings>"
} > "$out"
