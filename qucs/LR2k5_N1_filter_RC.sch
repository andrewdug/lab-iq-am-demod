<Qucs Schematic 25.2.99>
<Properties>
  <View=-50,-50,1800,1300,1,0,0>
  <Grid=10,10,1>
  <DataSet=LR2k5_N1_filter_RC.dat>
  <DataDisplay=LR2k5_N1_filter_RC.dpl>
  <OpenDisplay=0>
  <Script=LR2k5_N1_filter_RC.m>
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
  <Vac V1 1 250 270 18 -26 0 1 "1 V" 1 "1 kHz" 0 "0" 0 "0" 0 "0" 0 "0" 0>
  <R R1 1 330 240 -26 15 0 0 "1 kOhm" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "european" 0>
  <C C1 1 400 270 17 -26 0 1 "17.44 pF" 1 "" 0 "neutral" 0>
  <GND * 1 400 300 0 0 0 0>
  <.AC AC1 1 0 400 0 78 0 0 "log" 1 "100 kHz" 1 "1 GHz" 1 "801" 1 "no" 0>
  <Eqn Eqn1 1 300 450 -37 18 0 0 "KdB=dB(out.v)" 1 "tau=-diff(unwrap(angle(out.v)),2*pi*acfrequency)" 1 "yes" 0>
</Components>
<Wires>
  <360 240 400 240 "" 0 0 0 "">
  <250 240 300 240 "" 0 0 0 "">
  <250 300 400 300 "" 0 0 0 "">
  <400 240 400 240 "out" 430 210 0 "">
</Wires>
<Diagrams>
  <Rect 700 450 800 350 3 #c0c0c0 1 10 1 0 1 1 1 -100 20 10 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "K, дБ" "">
	<"KdB" #0000ff 2 3 0 0 0>
  </Rect>
  <Rect 700 950 800 350 3 #c0c0c0 1 10 1 0 1 1 1 0 1e-08 1e-07 1 -1 0.5 1 315 0 225 1 0 0 "Частота, Гц" "tau, с" "">
	<"tau" #ff0000 2 3 0 0 0>
  </Rect>
</Diagrams>
<Paintings>
  <Text 0 0 14 #00007f 0 "RC-ФНЧ бригады N = 1: R = 1 кОм, C = 17,44 пФ (tau = 17,44 нс), fc = 9,13 МГц">
  <Text 0 60 11 #7f0000 0 "Расчёт: K(3 МГц) = -0,45 дБ (0,95); K(35/40/45 МГц) = -12,0 / -13,1 / -14,0 дБ; K(860 МГц) = -39,5 дБ.\nГрупповая задержка на НЧ = RC = 17,4 нс.">
  <Text 700 50 14 #00007f 0 "АЧХ, дБ (сверху) и групповая задержка (снизу), 100 кГц - 1 ГГц">
</Paintings>
