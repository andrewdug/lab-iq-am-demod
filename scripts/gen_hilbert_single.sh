#!/bin/bash
# Qucs-S schematic with ONE FIR Hilbert filter, drawn after MicroCap ex5:
#   delay line of 2Dt blocks (Ideal/VDelay), coefficient gains on the taps (VCVS, T = 0),
#   two group adders (like X4, X12 in ex5; Ideal/VSum chains) and a final adder (like X7),
#   alignment delay X1 = t0 = (n0 - 1) Dt of the source signal (not part of the filter).
# Usage: gen_hilbert_single.sh <AC|TR> <R10|W10|W18> <out.sch>
#   R10 - order 10, rectangular window; W10 - order 10, Hamming; W18 - order 18, Hamming.
set -euo pipefail
mode=${1:?mode AC|TR}; filt=${2:?filter R10|W10|W18}; out=${3:?output .sch}; name=$(basename "$out" .sch)
case $mode in AC|TR) ;; *) echo "mode must be AC or TR" >&2; exit 1;; esac
case $filt in
  R10) L=10; win=0; title="порядок 10, без окна (прямоугольное окно)";;
  W10) L=10; win=1; title="порядок 10, окно Хэмминга";;
  W18) L=18; win=1; title="порядок 18, окно Хэмминга";;
  *) echo "filter must be R10, W10 or W18" >&2; exit 1;;
esac
nt=$((L / 2 + 1)); n0=$((L / 2 + 1)); d0=$((L / 2)); n1=$(((nt + 1) / 2))

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
comps=$tmp/comps; wires=$tmp/wires; diags=$tmp/diags
: > "$comps"; : > "$wires"; : > "$diags"
ROW=300; AY=$((ROW + 70)); YB=$((ROW + 130))     # gain row, adder row, adder link row

wire() {
  if [ -n "${5:-}" ]; then
    echo "  <$1 $2 $3 $4 \"$5\" $(( ($1 + $3) / 2 - 10 )) $(( $2 - 30 )) 0 \"\">" >> "$wires"
  else
    echo "  <$1 $2 $3 $4 \"\" 0 0 0 \"\">" >> "$wires"
  fi
}
gnd()   { echo "  <GND * 1 $1 $2 0 0 0 0>" >> "$comps"; }
text()  { echo "  <Text $1 $2 $3 #00007f 0 \"$4\">" >> "$diags"; }
delay() { echo "  <Lib $1 1 $2 $3 -20 24 0 0 \"Ideal\" 0 \"VDelay\" 0 \"$4\" 1>" >> "$comps"; }
adder() { echo "  <Lib $1 1 $2 $3 22 -36 0 0 \"Ideal\" 0 \"VSum\" 0>" >> "$comps"; }
gain()  { echo "  <VCVS $1 1 $2 $3 -26 34 0 0 \"$4\" 1 \"0\" 0>" >> "$comps"; gnd $(($2 - 30)) $(($3 + 30)); gnd $(($2 + 30)) $(($3 + 30)); }
rload() { echo "  <R $1 1 $2 $3 15 -26 0 1 \"1 MOhm\" 1 \"26.85\" 0 \"0.0\" 0 \"0.0\" 0 \"26.85\" 0 \"european\" 0>" >> "$comps"; gnd "$2" $(($3 + 30)); }

# --- source ---
if [ "$mode" = TR ]; then
  echo '  <AM_Mod V1 1 100 200 18 -26 0 0 "U0" 1 "f0" 1 "0" 1 "m" 1>' >> "$comps"
  echo '  <Vac Vm 1 -20 230 18 -26 0 1 "1" 1 "Fm" 1 "0" 0 "0" 0 "0" 0 "0" 0>' >> "$comps"
  wire -20 200 70 200; wire -20 260 100 260; wire 100 230 100 260; gnd 100 260
else
  echo '  <Vac V1 1 100 200 18 -26 0 1 "1 V" 1 "1 GHz" 0 "0" 0 "0" 0 "0" 0 "0" 0>' >> "$comps"
  gnd 100 230
fi
wire 100 170 100 70

# --- delay line: L/2 blocks of 2Dt, nodes x0, x2, ..., xL ---
wire 100 70 350 70 x0
for ((i = 1; i <= d0; i++)); do
  x=$((200 + 200 * i))
  delay "X$i" $x 70 T2
  wire $((x - 50)) 70 $((x - 30)) 70
  wire $((x + 30)) 70 $((x + 150)) 70 "x$((2 * i))"
done

# --- taps: bus down to the gain row; gain K_j = h(2j+1) ---
for ((j = 0; j < nt; j++)); do
  bx=$((350 + 200 * j)); x=$((400 + 200 * j)); ax=$((x + 100))
  wire $bx 70 $bx $((ROW - 30))
  gain "K$j" $x $ROW "c$j"
  wire $bx $((ROW - 30)) $((x - 30)) $((ROW - 30))
  wire $((x + 30)) $((ROW - 30)) $ax $((ROW - 30))
  if [ $j -eq 0 ] || [ $j -eq $n1 ]; then          # first tap of a group -> lower input of the next adder
    wire $ax $((ROW - 30)) $ax $YB
    wire $ax $YB $((ax + 200)) $YB
    wire $((ax + 200)) $YB $((ax + 200)) $((AY + 30))
  else                                             # adder A_j: upper input = K_j x(2j)
    adder "A$j" $ax $AY
    wire $ax $((ROW - 30)) $ax $((AY - 30))
    if [ $j -ne $((n1 - 1)) ] && [ $j -ne $((nt - 1)) ]; then   # chain inside the group
      wire $((ax + 30)) $AY $((ax + 30)) $YB
      wire $((ax + 30)) $YB $((ax + 200)) $YB
      wire $((ax + 200)) $YB $((ax + 200)) $((AY + 30))
    fi
  fi
done

# --- final adder (like X7): group 1 = taps 0..n1-1, group 2 = taps n1..nt-1 ---
g1=$((530 + 200 * (n1 - 1))); g2=$((530 + 200 * (nt - 1)))
fx=$((g2 + 100)); fy=$((ROW + 260))
adder Asum $fx $fy
wire $g1 $AY $g1 $((ROW + 200)); wire $g1 $((ROW + 200)) $fx $((ROW + 200)); wire $fx $((ROW + 200)) $fx $((fy - 30))
wire $g2 $AY $((fx - 40)) $AY; wire $((fx - 40)) $AY $((fx - 40)) $((fy + 60))
wire $((fx - 40)) $((fy + 60)) $fx $((fy + 60)); wire $fx $((fy + 60)) $fx $((fy + 30))
wire $((fx + 30)) $fy $((fx + 130)) $fy y
rload Ry $((fx + 130)) $((fy + 30))

# --- alignment delay X0 = t0 (like X1 in ex5), fed from the x0 bus ---
ry=$((ROW + 420))
wire 350 $((ROW - 30)) 350 $ry
delay X0 450 $ry T0
wire 350 $ry 420 $ry
wire 480 $ry 600 $ry xd
rload Rxd 600 $((ry + 30))

# --- coefficients ---
ceq=""; cnum=""
for ((j = 0; j < nt; j++)); do
  k=$((2 * j - d0))
  if [ $win = 1 ]; then expr="2/(pi*($k))*(0.54-0.46*cos(2*pi*$((2 * j))/$L))"; else expr="2/(pi*($k))"; fi
  ceq="$ceq \"c$j=$expr\" 1"
  cnum="$cnum$(awk -v k=$k -v n=$((2 * j)) -v L=$L -v w=$win 'BEGIN{pi=atan2(0,-1); c=2/(pi*k); if(w) c*=0.54-0.46*cos(2*pi*n/L); printf "%s%.4f", (n?" ":""), c}')"
done
cnum=$(echo "$cnum" | tr . ,)

cat >> "$comps" <<EOF
  <Eqn Eqn1 1 -300 -600 -37 18 0 0 "N=1" 1 "f0=(400+30*N)*1e6" 1 "Fm=3*N*1e6" 1 "U0=1" 1 "m=1" 1 "fs=4*f0" 1 "Dt=1/fs" 1 "T2=2*Dt" 1 "T0=$d0*Dt" 1 "yes" 0>
  <Eqn EqnC 1 300 -600 -37 18 0 0$ceq "yes" 0>
EOF

xr=$((fx + 400))
if [ "$mode" = TR ]; then
  cat >> "$comps" <<'EOF'
  <.TR TR1 1 -300 900 0 78 0 0 "lin" 1 "0" 1 "2e-6" 1 "41281" 1 "Trapezoidal" 0 "2" 0 "1 ps" 0 "1e-16" 0 "150" 0 "0.001" 0 "1 pA" 0 "1 uV" 0 "26.85" 0 "1e-3" 0 "1e-6" 0 "1" 0 "CroutLU" 0 "no" 0 "yes" 0 "10 ps" 0>
  <Eqn EqnE 1 900 -600 -37 18 0 0 "Uenv=U0*(1+m*sin(2*pi*Fm*(time-T0)))" 1 "Uz=sqrt(xd.Vt^2+y.Vt^2)" 1 "yes" 0>
EOF
  cat >> "$diags" <<EOF
  <Rect $xr -300 1000 300 3 #c0c0c0 1 00 1 0 1e-07 1e-06 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "t, с" "В" "">
	<"x0.Vt" #008000 0 3 0 0 0>
  </Rect>
  <Rect $xr 200 1000 300 3 #c0c0c0 1 00 0 1e-06 2e-09 1.01e-06 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "t, с" "В" "">
	<"x0.Vt" #008000 0 3 0 0 0>
	<"xd.Vt" #000000 1 3 0 0 0>
	<"y.Vt" #0000ff 1 3 0 0 0>
  </Rect>
  <Rect $xr 700 1000 300 3 #c0c0c0 1 00 1 0 1e-07 1e-06 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "t, с" "В" "">
	<"Uenv" #000000 2 3 0 1 0>
	<"Uz" #0000ff 0 3 0 0 0>
  </Rect>
  <Text $xr -380 14 #00007f 0 "Входной АМ-сигнал s(t) = x0: f0 = 430 МГц, Fм = 3 МГц, m = 1, U0 = 1 В">
  <Text $xr 120 14 #00007f 0 "s(t) (зелёный), s'(t) = s(t - t0) = xd (чёрный), сопряжённый по Гильберту s^(t) = y (синий)">
  <Text $xr 620 14 #00007f 0 "Огибающая sqrt(xd^2 + y^2) (синий) и модулирующее U0(1 + m sin 2pi Fм (t - t0)) (чёрный пунктир)">
EOF
else
  cat >> "$comps" <<'EOF'
  <.AC AC1 1 -300 900 0 78 0 0 "lin" 1 "1 MHz" 1 "860 MHz" 1 "860" 1 "no" 0>
  <Eqn EqnA 1 900 -600 -37 18 0 0 "K=mag(y.v)" 1 "Kd=mag(xd.v)" 1 "F=unwrap(angle(y.v))" 1 "Fd=unwrap(angle(xd.v))" 1 "dF=angle(y.v/xd.v)*180/pi" 1 "tg=-diff(F,acfrequency)/(2*pi)" 1 "K_f0=K[429]" 1 "dK_f0=(K[429]-Kd[429])*100" 1 "dF_f0=dF[429]" 1 "t0_ns=tg[429]*1e9" 1 "T0_ns=T0*1e9" 1 "yes" 0>
EOF
  cat >> "$diags" <<EOF
  <Rect $xr -300 1000 350 3 #c0c0c0 1 00 0 0 1e+08 8.6e+08 0 0 0.2 1.4 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "|K|" "">
	<"Kd" #000000 1 3 0 1 0>
	<"K" #0000ff 2 3 0 0 0>
  </Rect>
  <Rect $xr 250 1000 350 3 #c0c0c0 1 00 0 0 1e+08 8.6e+08 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "рад" "">
	<"Fd" #000000 1 3 0 1 0>
	<"F" #0000ff 2 3 0 0 0>
  </Rect>
  <Rect $xr 800 1000 350 3 #c0c0c0 1 00 0 0 1e+08 8.6e+08 0 -180 30 0 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "град" "">
	<"dF" #0000ff 2 3 0 0 0>
  </Rect>
  <Tab $xr 1300 500 200 3 #c0c0c0 1 00 1 0 1 1 1 0 1 1 1 0 1 1 315 0 225 0 0 0 "" "" "">
	<"K_f0" #0000ff 0 4 0 0 0>
	<"dK_f0" #0000ff 0 3 0 0 0>
	<"dF_f0" #0000ff 0 3 0 0 0>
	<"t0_ns" #0000ff 0 4 0 0 0>
	<"T0_ns" #0000ff 0 4 0 0 0>
  </Tab>
  <Text $xr -380 14 #00007f 0 "АЧХ фильтра (синий) и элемента задержки X0 = t0 (чёрный пунктир); f0 = fs/4 = 430 МГц">
  <Text $xr 170 14 #00007f 0 "ФЧХ (развёрнутая), рад: фильтр (синий), элемент задержки (чёрный пунктир)">
  <Text $xr 720 14 #00007f 0 "Разность фаз выхода фильтра и элемента задержки X0, град">
  <Text $xr 1220 14 #00007f 0 "На f0 = 430 МГц: |K|, отклонение от задержки (%), разность фаз (град), задержка по ФЧХ и теоретическая (нс)">
EOF
fi

text -300 -900 14 "Цифровой КИХ-фильтр Гильберта, $title (N = $((L + 1)), n0 = $n0), по образцу MicroCap ex5. N = 1: fs = 4 f0 = 1720 МГц, Dt = 0,5814 нс.\nКоэффициенты h(n) = 2/(pi (n - n0)) при нечётном n - n0$( [ $win = 1 ] && echo ', умноженные на окно Хэмминга w(n) = 0,54 - 0,46 cos(2 pi (n-1)/(N-1))' ); ненулевые: $cnum.\nЗадержка фильтра t0 = (n0 - 1) Dt = ${d0}Dt."
text -300 -110 14 "ЛИНИЯ ЗАДЕРЖКИ: X1...X$d0 - блоки VDelay (библиотека Ideal), Tdel = 2Dt; узлы x0, x2, ..., x$L = s(t - k Dt)."
text -300 $((ROW - 60)) 14 "УСИЛИТЕЛИ K0...K$((nt - 1)):\nИНУН без задержки, K = c0...c$((nt - 1))\n(нулевые коэффициенты пропущены)"
text -300 $((ROW + 120)) 14 "СУММАТОРЫ (VSum, Ideal):\nгруппа 1 - отводы 0...$((n1 - 1)), группа 2 - отводы $n1...$((nt - 1)),\nитоговый Asum -> y = s^(t)"
text -300 $((ry - 60)) 14 "X0 - ВЫРАВНИВАЮЩАЯ ЗАДЕРЖКА t0 = ${d0}Dt\n(как X1 в ex5, в фильтр не входит) -> xd = s(t - t0)"

{
cat <<EOF
<Qucs Schematic 25.2.99>
<Properties>
  <View=-400,-950,$((xr + 1100)),1500,0.5,0,0>
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
