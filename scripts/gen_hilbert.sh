#!/bin/bash
# Generates Qucs-S schematics of a 10th-order FIR Hilbert filter (MicroCap ex5 analogue)
# Usage: gen_hilbert.sh <mode: AC|TR> <out.sch>
mode=$1; out=$2; name=$(basename "$out" .sch)

lbl() { echo "  <$1 $2 $1 $2 \"$3\" $(( $1 + 10 )) $(( $2 - 30 )) 0 \"\">" >> wires; }
gnd() { echo "  <GND * 1 $1 $2 0 0 0 0>" >> comps; }
# vcvs NAME X Y G T IN OUTP OUTN   (IN/OUTP/OUTN: node label or "gnd"; in- always gnd)
vcvs() {
  local n=$1 x=$2 y=$3 g=$4 t=$5 in=$6 op=$7 on=$8
  echo "  <VCVS $n 1 $x $y -26 34 0 0 \"$g\" 1 \"$t\" 1>" >> comps
  lbl $((x-30)) $((y-30)) "$in"
  gnd $((x-30)) $((y+30))
  [ "$op" = gnd ] && gnd $((x+30)) $((y-30)) || lbl $((x+30)) $((y-30)) "$op"
  [ "$on" = gnd ] && gnd $((x+30)) $((y+30)) || lbl $((x+30)) $((y+30)) "$on"
}
rload() { echo "  <R $1 1 $2 $3 15 -26 0 1 \"1 MOhm\" 1 \"26.85\" 0 \"0.0\" 0 \"0.0\" 0 \"26.85\" 0 \"european\" 0>" >> comps
  lbl $2 $(($3-30)) "$4"; gnd $2 $(($3+30)); }

: > comps; : > wires; : > diags

# --- source ---
if [ "$mode" = TR ]; then
  echo '  <AM_Mod V1 1 100 200 18 -26 0 0 "U0" 1 "f0" 1 "0" 1 "m" 1>' >> comps
  echo '  <Vac Vm 1 -20 230 18 -26 0 1 "1" 1 "Fm" 1 "0" 0 "0" 0 "0" 0 "0" 0>' >> comps
  echo '  <-20 200 70 200 "" 0 0 0 "">' >> wires
  echo '  <-20 260 100 260 "" 0 0 0 "">' >> wires
  echo '  <100 230 100 260 "" 0 0 0 "">' >> wires
  gnd 100 260
else
  echo '  <Vac V1 1 100 200 18 -26 0 1 "1 V" 1 "1 GHz" 0 "0" 0 "0" 0 "0" 0 "0" 0>' >> comps
  gnd 100 230
fi
lbl 100 170 x0

# --- delay line X3..X14: five 2*Dt sections, plus 5*Dt alignment delay X1 ---
for i in 1 2 3 4 5; do vcvs D$i $((200+200*i)) 100 1 T2 x$((2*i-2)) x$((2*i)) gnd; done
vcvs D0 400 300 1 T5 x0 xd gnd

# --- taps: rectangular (H) and Hamming-weighted (W), outputs stacked in series = summator ---
for k in 0 1 2 3 4 5; do
  x=$((400+200*k)); on=gnd; [ $k -gt 0 ] && on=s$k; op=s$((k+1)); [ $k = 5 ] && op=y
  vcvs H$k $x 500 c$k 0 x$((2*k)) $op $on
  on=gnd; [ $k -gt 0 ] && on=sw$k; op=sw$((k+1)); [ $k = 5 ] && op=yw
  vcvs W$k $x 700 cw$k 0 x$((2*k)) $op $on
done
rload Ry 1700 500 y; rload Ryw 1800 500 yw; rload Rxd 1900 500 xd

cat >> comps <<'EOF'
  <Eqn Eqn1 1 -300 -400 -37 18 0 0 "N=1" 1 "f0=(400+30*N)*1e6" 1 "Fm=3*N*1e6" 1 "U0=1" 1 "m=1" 1 "fs=4*f0" 1 "Dt=1/fs" 1 "T2=2*Dt" 1 "T5=5*Dt" 1 "yes" 0>
  <Eqn EqnH 1 200 -400 -37 18 0 0 "c0=-2/(5*pi)" 1 "c1=-2/(3*pi)" 1 "c2=-2/pi" 1 "c3=2/pi" 1 "c4=2/(3*pi)" 1 "c5=2/(5*pi)" 1 "w5=0.54-0.46*cos(2*pi*0/10)" 1 "w3=0.54-0.46*cos(2*pi*2/10)" 1 "w1=0.54-0.46*cos(2*pi*4/10)" 1 "cw0=c0*w5" 1 "cw1=c1*w3" 1 "cw2=c2*w1" 1 "cw3=c3*w1" 1 "cw4=c4*w3" 1 "cw5=c5*w5" 1 "yes" 0>
EOF

if [ "$mode" = TR ]; then
  cat >> comps <<'EOF'
  <.TR TR1 1 -300 400 0 78 0 0 "lin" 1 "0" 1 "2e-6" 1 "40001" 1 "Trapezoidal" 0 "2" 0 "1 ps" 0 "1e-16" 0 "150" 0 "0.001" 0 "1 pA" 0 "1 uV" 0 "26.85" 0 "1e-3" 0 "1e-6" 0 "1" 0 "CroutLU" 0 "no" 0 "yes" 0 "10 ps" 0>
  <Eqn EqnE 1 700 -400 -37 18 0 0 "Uenv=U0*(1+m*sin(2*pi*Fm*(time-T5)))" 1 "Ur=sqrt(xd.Vt^2+y.Vt^2)" 1 "Uw=sqrt(xd.Vt^2+yw.Vt^2)" 1 "yes" 0>
EOF
  cat >> diags <<'EOF'
  <Rect 2200 -300 1000 300 3 #c0c0c0 1 00 1 0 1e-07 1e-06 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "t, с" "В" "">
	<"x0.Vt" #0000ff 0 3 0 0 0>
  </Rect>
  <Rect 2200 200 1000 300 3 #c0c0c0 1 00 0 1e-06 2e-09 1.01e-06 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "t, с" "В" "">
	<"xd.Vt" #000000 1 3 0 0 0>
	<"y.Vt" #ff0000 1 3 0 0 0>
	<"yw.Vt" #0000ff 1 3 0 1 0>
  </Rect>
  <Rect 2200 700 1000 300 3 #c0c0c0 1 00 1 0 1e-07 1e-06 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "t, с" "В" "">
	<"Uenv" #000000 2 3 0 0 0>
	<"Ur" #ff0000 0 3 0 0 0>
  </Rect>
  <Rect 2200 1200 1000 300 3 #c0c0c0 1 00 1 0 1e-07 1e-06 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "t, с" "В" "">
	<"Uenv" #000000 2 3 0 0 0>
	<"Uw" #0000ff 0 3 0 0 0>
  </Rect>
  <Text 2200 -380 14 #00007f 0 "Входной АМ-сигнал x0: f0 = 430 МГц, Fм = 3 МГц, m = 1, U0 = 1 В">
  <Text 2200 120 14 #00007f 0 "Выход элемента задержки xd (чёрный) и фильтра Гильберта: y - без окна (красный), yw - окно Хэмминга (синий пунктир)">
  <Text 2200 620 14 #00007f 0 "Огибающая sqrt(xd^2 + y^2) без окна (красный) и модулирующее U0(1+m sin 2pi Fм (t-5Dt)) (чёрный)">
  <Text 2200 1120 14 #00007f 0 "Огибающая sqrt(xd^2 + yw^2) с окном Хэмминга (синий) и модулирующее (чёрный)">
EOF
else
  cat >> comps <<'EOF'
  <.AC AC1 1 -300 400 0 78 0 0 "lin" 1 "1 MHz" 1 "860 MHz" 1 "860" 1 "no" 0>
  <Eqn EqnA 1 700 -400 -37 18 0 0 "Kh=mag(y.v)" 1 "Kw=mag(yw.v)" 1 "Kd=mag(xd.v)" 1 "Fh=unwrap(angle(y.v))" 1 "Fw=unwrap(angle(yw.v))" 1 "Fd=unwrap(angle(xd.v))" 1 "dFh=angle(y.v/xd.v)*180/pi" 1 "dFw=angle(yw.v/xd.v)*180/pi" 1 "Kh_f0=Kh[429]" 1 "Kw_f0=Kw[429]" 1 "dKh_f0=(Kh[429]-Kd[429])*100" 1 "dKw_f0=(Kw[429]-Kd[429])*100" 1 "dFh_f0=dFh[429]" 1 "dFw_f0=dFw[429]" 1 "yes" 0>
EOF
  cat >> diags <<'EOF'
  <Rect 2200 -300 1000 350 3 #c0c0c0 1 00 0 0 1e+08 8.6e+08 0 0 0.2 1.4 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "|K|" "">
	<"Kd" #000000 1 3 0 1 0>
	<"Kh" #ff0000 2 3 0 0 0>
	<"Kw" #0000ff 2 3 0 0 0>
  </Rect>
  <Rect 2200 250 1000 350 3 #c0c0c0 1 00 0 0 1e+08 8.6e+08 1 -1 0.5 1 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "рад" "">
	<"Fd" #000000 1 3 0 1 0>
	<"Fh" #ff0000 2 3 0 0 0>
	<"Fw" #0000ff 1 3 0 1 0>
  </Rect>
  <Rect 2200 800 1000 350 3 #c0c0c0 1 00 0 0 1e+08 8.6e+08 0 -180 30 0 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "град" "">
	<"dFh" #ff0000 2 3 0 0 0>
	<"dFw" #0000ff 1 3 0 1 0>
  </Rect>
  <Tab 2200 1300 600 200 3 #c0c0c0 1 00 1 0 1 1 1 0 1 1 1 0 1 1 315 0 225 0 0 0 "" "" "">
	<"Kh_f0" #0000ff 0 4 0 0 0>
	<"Kw_f0" #0000ff 0 4 0 0 0>
	<"dKh_f0" #0000ff 0 3 0 0 0>
	<"dKw_f0" #0000ff 0 3 0 0 0>
	<"dFh_f0" #0000ff 0 3 0 0 0>
	<"dFw_f0" #0000ff 0 3 0 0 0>
  </Tab>
  <Text 2200 -380 14 #00007f 0 "АЧХ: элемент задержки 5Dt (чёрный пунктир), фильтр Гильберта без окна (красный) и с окном Хэмминга (синий); f0 = fs/4 = 430 МГц">
  <Text 2200 170 14 #00007f 0 "ФЧХ (развёрнутая), рад: задержка (чёрный), Гильберт без окна (красный), с окном (синий)">
  <Text 2200 720 14 #00007f 0 "Разность фаз выходов фильтра Гильберта и элемента задержки, град">
  <Text 2200 1220 14 #00007f 0 "На частоте сигнала f0 = 430 МГц: |K|, отклонение от задержки (%), разность фаз (град)">
EOF
fi
cat >> diags <<'EOF'
  <Text -300 -600 14 #00007f 0 "Цифровой КИХ-фильтр Гильберта 10-го порядка (аналог модели MicroCap ex5), N = 1:\nfs = 4 f0 = 1720 МГц, Dt = 1/fs = 0,5814 нс; линия задержки 5 x 2Dt, элемент выравнивания 5Dt.\nh[k] = 2/(pi k) для нечётных k = -5..5 (чётные = 0): -0,1273 -0,2122 -0,6366 0 0,6366 0,2122 0,1273.\nОкно Хэмминга w[n] = 0,54 - 0,46 cos(2 pi n/(N-1)), N = 11: 0,08 0,398 0,912 1 0,912 0,398 0,08.\nСумматор - выходы ИНУН H0..H5 (W0..W5), включённые последовательно.">
EOF

{
cat <<EOF
<Qucs Schematic 25.2.99>
<Properties>
  <View=-400,-700,3300,1600,0.5,0,0>
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
cat comps
echo "</Components>"; echo "<Wires>"; cat wires; echo "</Wires>"
echo "<Diagrams>"; grep -v '^  <Text' diags; echo "</Diagrams>"; echo "<Paintings>"; grep '^  <Text' diags; echo "</Paintings>"
} > "$out"
rm -f comps wires diags
