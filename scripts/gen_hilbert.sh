#!/bin/bash
# Generates Qucs-S schematics of FIR Hilbert filters (MicroCap ex5 analogue):
#   H  - order 10 (N = 11, n0 = 6), rectangular window      -> node y,   reference xd   (5 Dt)
#   W  - order 10, Hamming window                            -> node yw,  reference xd   (5 Dt)
#   G  - order 18 (N = 19, n0 = 10), Hamming window (p. 4)   -> node y18, reference xd18 (9 Dt)
# Usage: gen_hilbert.sh <mode: AC|TR> <out.sch>
set -euo pipefail
mode=${1:?mode AC|TR}; out=${2:?output .sch}; name=$(basename "$out" .sch)
case $mode in AC|TR) ;; *) echo "mode must be AC or TR" >&2; exit 1;; esac

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
comps=$tmp/comps; wires=$tmp/wires; diags=$tmp/diags
: > "$comps"; : > "$wires"; : > "$diags"

lbl() { echo "  <$1 $2 $1 $2 \"$3\" $(( $1 + 10 )) $(( $2 - 30 )) 0 \"\">" >> "$wires"; }
gnd() { echo "  <GND * 1 $1 $2 0 0 0 0>" >> "$comps"; }
# vcvs NAME X Y G T IN OUTP OUTN   (IN/OUTP/OUTN: node label or "gnd"; in- always gnd)
vcvs() {
  local n=$1 x=$2 y=$3 g=$4 t=$5 in=$6 op=$7 on=$8
  echo "  <VCVS $n 1 $x $y -26 34 0 0 \"$g\" 1 \"$t\" 1>" >> "$comps"
  lbl $((x-30)) $((y-30)) "$in"
  gnd $((x-30)) $((y+30))
  if [ "$op" = gnd ]; then gnd $((x+30)) $((y-30)); else lbl $((x+30)) $((y-30)) "$op"; fi
  if [ "$on" = gnd ]; then gnd $((x+30)) $((y+30)); else lbl $((x+30)) $((y+30)) "$on"; fi
}
rload() { echo "  <R $1 1 $2 $3 15 -26 0 1 \"1 MOhm\" 1 \"26.85\" 0 \"0.0\" 0 \"0.0\" 0 \"26.85\" 0 \"european\" 0>" >> "$comps"
  lbl "$2" $(($3-30)) "$4"; gnd "$2" $(($3+30)); }
# chain PREFIX ROW ORDER COEF NODE OUT: taps x0, x2, ..., x<ORDER>; VCVS outputs in series = summator
chain() {
  local p=$1 row=$2 order=$3 cf=$4 nd=$5 outn=$6 nt=$(( $3 / 2 + 1 )) j op on
  for ((j = 0; j < nt; j++)); do
    on=gnd; [ $j -gt 0 ] && on=$nd$j
    op=$nd$((j + 1)); [ $j -eq $((nt - 1)) ] && op=$outn
    vcvs "$p$j" $((400 + 200 * j)) "$row" "$cf$j" 0 "x$((2 * j))" "$op" "$on"
  done
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
lbl 100 170 x0

# --- delay line D1..D9: nine 2*Dt sections (x0 ... x18); alignment delays 5*Dt (xd) and 9*Dt (xd18) ---
for i in 1 2 3 4 5 6 7 8 9; do vcvs D$i $((200 + 200 * i)) 100 1 T2 x$((2 * i - 2)) x$((2 * i)) gnd; done
vcvs D0 400 300 1 T5 x0 xd gnd
vcvs D18 800 300 1 T9 x0 xd18 gnd

# --- filters ---
chain H 500 10 c s y
chain W 700 10 cw sw yw
chain G 900 18 g sg y18
rload Ry 2400 500 y; rload Ryw 2500 500 yw; rload Ry18 2600 500 y18
rload Rxd 2700 500 xd; rload Rxd18 2800 500 xd18

# order-18 Hamming coefficients: g_j = 2/(pi k) * (0.54 - 0.46 cos(2 pi n / 18)), n = 2j, k = n - 9
geq=""
for j in 0 1 2 3 4 5 6 7 8 9; do
  geq="$geq \"g$j=2/(pi*($((2 * j - 9))))*(0.54-0.46*cos(2*pi*$((2 * j))/18))\" 1"
done
gnum=$(awk 'BEGIN{pi=atan2(0,-1); for(j=0;j<10;j++){n=2*j; k=n-9; printf "%s%.4f", (j?" ":""), 2/(pi*k)*(0.54-0.46*cos(2*pi*n/18))}}' | tr . ,)

cat >> "$comps" <<EOF
  <Eqn Eqn1 1 -300 -400 -37 18 0 0 "N=1" 1 "f0=(400+30*N)*1e6" 1 "Fm=3*N*1e6" 1 "U0=1" 1 "m=1" 1 "fs=4*f0" 1 "Dt=1/fs" 1 "T2=2*Dt" 1 "T5=5*Dt" 1 "T9=9*Dt" 1 "yes" 0>
  <Eqn EqnH 1 200 -400 -37 18 0 0 "c0=-2/(5*pi)" 1 "c1=-2/(3*pi)" 1 "c2=-2/pi" 1 "c3=2/pi" 1 "c4=2/(3*pi)" 1 "c5=2/(5*pi)" 1 "w5=0.54-0.46*cos(2*pi*0/10)" 1 "w3=0.54-0.46*cos(2*pi*2/10)" 1 "w1=0.54-0.46*cos(2*pi*4/10)" 1 "cw0=c0*w5" 1 "cw1=c1*w3" 1 "cw2=c2*w1" 1 "cw3=c3*w1" 1 "cw4=c4*w3" 1 "cw5=c5*w5" 1 "yes" 0>
  <Eqn EqnG 1 700 -400 -37 18 0 0$geq "yes" 0>
EOF

if [ "$mode" = TR ]; then
  # 41281 points: step Dt/12, so the samples t = n*Dt lie exactly on the output grid
  cat >> "$comps" <<'EOF'
  <.TR TR1 1 -300 400 0 78 0 0 "lin" 1 "0" 1 "2e-6" 1 "41281" 1 "Trapezoidal" 0 "2" 0 "1 ps" 0 "1e-16" 0 "150" 0 "0.001" 0 "1 pA" 0 "1 uV" 0 "26.85" 0 "1e-3" 0 "1e-6" 0 "1" 0 "CroutLU" 0 "no" 0 "yes" 0 "10 ps" 0>
  <Eqn EqnE 1 1200 -400 -37 18 0 0 "Uenv=U0*(1+m*sin(2*pi*Fm*(time-T5)))" 1 "Uenv18=U0*(1+m*sin(2*pi*Fm*(time-T9)))" 1 "Ur=sqrt(xd.Vt^2+y.Vt^2)" 1 "Uw=sqrt(xd.Vt^2+yw.Vt^2)" 1 "Ug=sqrt(xd18.Vt^2+y18.Vt^2)" 1 "yes" 0>
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
  <Text 3000 -380 14 #00007f 0 "s(t) = x0: входной АМ-сигнал, f0 = 430 МГц, Fм = 3 МГц, m = 1, U0 = 1 В">
  <Text 3000 120 14 #00007f 0 "Порядок 10: s(t) (зелёный), s'(t) = s(t - 5Dt) = xd (чёрный), сопряжённый y - без окна (красный), yw - окно Хэмминга (синий пунктир)">
  <Text 3000 620 14 #00007f 0 "Порядок 18, окно Хэмминга: s(t) (зелёный), s'(t) = s(t - 9Dt) = xd18 (чёрный), сопряжённый y18 (синий)">
  <Text 3000 1120 14 #00007f 0 "Порядок 10: огибающая sqrt(xd^2 + y^2) без окна (красный), sqrt(xd^2 + yw^2) с окном (синий), модулирующее U0(1 + m sin 2pi Fм (t - 5Dt)) (чёрный)">
  <Text 3000 1620 14 #00007f 0 "Порядок 18: огибающая sqrt(xd18^2 + y18^2) (синий) и модулирующее U0(1 + m sin 2pi Fм (t - 9Dt)) (чёрный)">
EOF
else
  cat >> "$comps" <<'EOF'
  <.AC AC1 1 -300 400 0 78 0 0 "lin" 1 "1 MHz" 1 "860 MHz" 1 "860" 1 "no" 0>
  <Eqn EqnA 1 1200 -400 -37 18 0 0 "Kh=mag(y.v)" 1 "Kw=mag(yw.v)" 1 "Kg=mag(y18.v)" 1 "Kd=mag(xd.v)" 1 "Fh=unwrap(angle(y.v))" 1 "Fw=unwrap(angle(yw.v))" 1 "Fg=unwrap(angle(y18.v))" 1 "Fd=unwrap(angle(xd.v))" 1 "dFh=angle(y.v/xd.v)*180/pi" 1 "dFw=angle(yw.v/xd.v)*180/pi" 1 "dFg=angle(y18.v/xd18.v)*180/pi" 1 "tgw=-diff(Fw,acfrequency)/(2*pi)" 1 "tgg=-diff(Fg,acfrequency)/(2*pi)" 1 "yes" 0>
  <Eqn EqnF0 1 1700 -400 -37 18 0 0 "Kh_f0=Kh[429]" 1 "Kw_f0=Kw[429]" 1 "Kg_f0=Kg[429]" 1 "dKh_f0=(Kh[429]-Kd[429])*100" 1 "dKw_f0=(Kw[429]-Kd[429])*100" 1 "dKg_f0=(Kg[429]-Kd[429])*100" 1 "dFh_f0=dFh[429]" 1 "dFw_f0=dFw[429]" 1 "dFg_f0=dFg[429]" 1 "t0w_ns=tgw[429]*1e9" 1 "t0g_ns=tgg[429]*1e9" 1 "T5_ns=T5*1e9" 1 "T9_ns=T9*1e9" 1 "yes" 0>
EOF
  cat >> "$diags" <<'EOF'
  <Rect 3000 -300 1000 350 3 #c0c0c0 1 00 0 0 1e+08 8.6e+08 0 0 0.2 1.4 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "|K|" "">
	<"Kd" #000000 1 3 0 1 0>
	<"Kh" #ff0000 2 3 0 0 0>
	<"Kw" #0000ff 2 3 0 0 0>
	<"Kg" #008000 2 3 0 0 0>
  </Rect>
  <Rect 3000 250 1000 350 3 #c0c0c0 1 00 0 0 1e+08 8.6e+08 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "рад" "">
	<"Fd" #000000 1 3 0 1 0>
	<"Fh" #ff0000 2 3 0 0 0>
	<"Fw" #0000ff 1 3 0 1 0>
	<"Fg" #008000 1 3 0 0 0>
  </Rect>
  <Rect 3000 800 1000 350 3 #c0c0c0 1 00 0 0 1e+08 8.6e+08 0 -180 30 0 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "град" "">
	<"dFh" #ff0000 2 3 0 0 0>
	<"dFw" #0000ff 1 3 0 1 0>
	<"dFg" #008000 1 3 0 0 0>
  </Rect>
  <Rect 3000 1350 1000 350 3 #c0c0c0 1 00 0 1e+08 1e+08 7.6e+08 0 0 1e-09 6e-09 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "с" "">
	<"tgw" #0000ff 1 3 0 0 0>
	<"tgg" #008000 1 3 0 0 0>
  </Rect>
  <Tab 3000 1850 700 300 3 #c0c0c0 1 00 1 0 1 1 1 0 1 1 1 0 1 1 315 0 225 0 0 0 "" "" "">
	<"Kh_f0" #0000ff 0 4 0 0 0>
	<"Kw_f0" #0000ff 0 4 0 0 0>
	<"Kg_f0" #0000ff 0 4 0 0 0>
	<"dKh_f0" #0000ff 0 3 0 0 0>
	<"dKw_f0" #0000ff 0 3 0 0 0>
	<"dKg_f0" #0000ff 0 3 0 0 0>
	<"dFh_f0" #0000ff 0 3 0 0 0>
	<"dFw_f0" #0000ff 0 3 0 0 0>
	<"dFg_f0" #0000ff 0 3 0 0 0>
	<"t0w_ns" #0000ff 0 4 0 0 0>
	<"T5_ns" #0000ff 0 4 0 0 0>
	<"t0g_ns" #0000ff 0 4 0 0 0>
	<"T9_ns" #0000ff 0 4 0 0 0>
  </Tab>
  <Text 3000 -380 14 #00007f 0 "АЧХ: задержка (чёрный пунктир), порядок 10 без окна (красный), порядок 10 с окном Хэмминга (синий), порядок 18 с окном Хэмминга (зелёный); f0 = fs/4 = 430 МГц">
  <Text 3000 170 14 #00007f 0 "ФЧХ (развёрнутая), рад: задержка 5Dt (чёрный), порядок 10 без окна (красный), с окном (синий), порядок 18 с окном (зелёный)">
  <Text 3000 720 14 #00007f 0 "Разность фаз выхода фильтра и своего элемента задержки (5Dt для порядка 10, 9Dt для порядка 18), град">
  <Text 3000 1270 14 #00007f 0 "Групповая задержка tгр = -dФ/dw: порядок 10 (синий, теория 5Dt = 2,907 нс), порядок 18 (зелёный, теория 9Dt = 5,233 нс)">
  <Text 3000 1770 14 #00007f 0 "На частоте f0 = 430 МГц: |K|, отклонение |K| от 1 (%), разность фаз (град), задержка по ФЧХ и теоретическая (нс)">
EOF
fi
cat >> "$diags" <<EOF
  <Text -300 -650 14 #00007f 0 "Цифровые КИХ-фильтры Гильберта (аналог модели MicroCap ex5), N = 1:\nfs = 4 f0 = 1720 МГц (f0 - в середине области равномерности АЧХ), Dt = 1/fs = 0,5814 нс; линия задержки 9 x 2Dt.\nПорядок 10 (11 коэфф., n0 = 6): h[k] = 2/(pi k), k = -5..5 нечётные: -0,1273 -0,2122 -0,6366 0 0,6366 0,2122 0,1273; задержка t0 = 5Dt = 2,907 нс.\nОкно Хэмминга w[n] = 0,54 - 0,46 cos(2 pi n/(N-1)), N = 11: 0,08 0,398 0,912 1 0,912 0,398 0,08.\nПорядок 18 (19 коэфф., n0 = 10) с окном Хэмминга, k = -9..9 нечётные: $gnum; t0 = 9Dt = 5,233 нс.\nСумматор - выходы ИНУН, включённые последовательно.">
EOF

{
cat <<EOF
<Qucs Schematic 25.2.99>
<Properties>
  <View=-400,-700,4100,2200,0.5,0,0>
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
