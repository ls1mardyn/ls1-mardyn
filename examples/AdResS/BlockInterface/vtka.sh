#!/usr/bin/zsh
VTKA=../../../cmake-build-debug/tools/VTKAnalysis
InputDir=../../src/

#(cd $VTKA && ./vtk-analysis -n 400 -w 20.0 --bbox0 200.0 --bbox1 200.0 --bbox2 200.0 -o ./vtka_wEUC.txt "$InputDir/AdResS_wEUC_node0_100000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 400 -w 20.0 --bbox0 200.0 --bbox1 200.0 --bbox2 200.0 -o ./vtka_wMAN.txt "$InputDir/AdResS_wMAN_node0_100000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 400 -w 20.0 --bbox0 200.0 --bbox1 200.0 --bbox2 200.0 -o ./vtka_wCOM.txt "$InputDir/AdResS_wCOM_node0_100000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 400 -w 20.0 --bbox0 200.0 --bbox1 200.0 --bbox2 200.0 -o ./vtka_wNEA.txt "$InputDir/AdResS_wNEA_node0_100000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 400 -w 20.0 --bbox0 200.0 --bbox1 200.0 --bbox2 200.0 -o ./vtka_wOFF.txt "$InputDir/AdResS_wOFF_node0_100000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 400 -w 20.0 --bbox0 200.0 --bbox1 200.0 --bbox2 200.0 -o ./vtka_wREF.txt "$InputDir/AdResS_wREF_NOAdResS_node0_100000.vtu" )
#exit 0

## 1k
#echo "[Progress] 0/6"
#(cd $VTKA && ./vtk-analysis -n 300 -w 10.0 --bbox0 300.0 --bbox1 300.0 --bbox2 300.0 -o ./vtka_CH4_1k.txt "$InputDir/BlockInterface_1k_CH4_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 300 -w 10.0 --bbox0 300.0 --bbox1 300.0 --bbox2 300.0 -o ./vtka_CH4_1k_REF.txt "$InputDir/BlockInterface_1k_CH4_REF_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 180 -w 10.0 --bbox0 180.0 --bbox1 180.0 --bbox2 180.0 -o ./vtka_C6H12_1k.txt "$InputDir/BlockInterface_1k_C6H12_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 180 -w 10.0 --bbox0 180.0 --bbox1 180.0 --bbox2 180.0 -o ./vtka_C6H12_1k_REF.txt "$InputDir/BlockInterface_1k_C6H12_REF_node0_10000.vtu" )
## 4k
#echo "[Progress] 1/6"
#(cd $VTKA && ./vtk-analysis -n 470 -w 10.0 --bbox0 470.0 --bbox1 470.0 --bbox2 460.0 -o ./vtka_CH4_4k.txt "$InputDir/BlockInterface_4k_CH4_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 470 -w 10.0 --bbox0 470.0 --bbox1 470.0 --bbox2 460.0 -o ./vtka_CH4_4k_REF.txt "$InputDir/BlockInterface_4k_CH4_REF_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 280 -w 10.0 --bbox0 280.0 --bbox1 280.0 --bbox2 280.0 -o ./vtka_C6H12_4k.txt "$InputDir/BlockInterface_4k_C6H12_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 280 -w 10.0 --bbox0 280.0 --bbox1 280.0 --bbox2 280.0 -o ./vtka_C6H12_4k_REF.txt "$InputDir/BlockInterface_4k_C6H12_REF_node0_10000.vtu" )
## 16k
#echo "[Progress] 2/6"
#(cd $VTKA && ./vtk-analysis -n 740 -w 10.0 --bbox0 740.0 --bbox1 740.0 --bbox2 740.0 -o ./vtka_CH4_16k.txt "$InputDir/BlockInterface_16k_CH4_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 740 -w 10.0 --bbox0 740.0 --bbox1 740.0 --bbox2 740.0 -o ./vtka_CH4_16k_REF.txt "$InputDir/BlockInterface_16k_CH4_REF_node0_10000.vtu" )
(cd $VTKA && ./vtk-analysis -n 450 -w 10.0 --bbox0 450.0 --bbox1 450.0 --bbox2 450.0 -o ./vtka_fth_C6H12_16k.txt "$InputDir/BlockInterface_16k_C6H12_node0_100000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 450 -w 10.0 --bbox0 450.0 --bbox1 450.0 --bbox2 450.0 -o ./vtka_C6H12_16k_REF.txt "$InputDir/BlockInterface_16k_C6H12_REF_node0_10000.vtu" )
## 64k
#echo "[Progress] 3/6"
#(cd $VTKA && ./vtk-analysis -n 1180 -w 10.0 --bbox0 1180.0 --bbox1 1180.0 --bbox2 1180.0 -o ./vtka_CH4_64k.txt "$InputDir/BlockInterface_64k_CH4_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 1180 -w 10.0 --bbox0 1180.0 --bbox1 1180.0 --bbox2 1180.0 -o ./vtka_CH4_64k_REF.txt "$InputDir/BlockInterface_64k_CH4_REF_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 720 -w 10.0 --bbox0 720.0 --bbox1 720.0 --bbox2 720.0 -o ./vtka_C6H12_64k.txt "$InputDir/BlockInterface_64k_C6H12_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 720 -w 10.0 --bbox0 720.0 --bbox1 720.0 --bbox2 720.0 -o ./vtka_C6H12_64k_REF.txt "$InputDir/BlockInterface_64k_C6H12_REF_node0_10000.vtu" )
## 256k
#echo "[Progress] 4/6"
#(cd $VTKA && ./vtk-analysis -n 1870 -w 10.0 --bbox0 1870.0 --bbox1 1870.0 --bbox2 1870.0 -o ./vtka_CH4_256k.txt "$InputDir/BlockInterface_256k_CH4_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 1870 -w 10.0 --bbox0 1870.0 --bbox1 1870.0 --bbox2 1870.0 -o ./vtka_CH4_256k_REF.txt "$InputDir/BlockInterface_256k_CH4_REF_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 1130 -w 10.0 --bbox0 1130.0 --bbox1 1130.0 --bbox2 1130.0 -o ./vtka_C6H12_256k.txt "$InputDir/BlockInterface_256k_C6H12_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 1130 -w 10.0 --bbox0 1130.0 --bbox1 1130.0 --bbox2 1130.0 -o ./vtka_C6H12_256k_REF.txt "$InputDir/BlockInterface_256k_C6H12_REF_node0_10000.vtu" )
## 1024k
#echo "[Progress] 5/6"
#(cd $VTKA && ./vtk-analysis -n 2970 -w 10.0 --bbox0 2970.0 --bbox1 2970.0 --bbox2 2970.0 -o ./vtka_CH4_1024k.txt "$InputDir/BlockInterface_1024k_CH4_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 2970 -w 10.0 --bbox0 2970.0 --bbox1 2970.0 --bbox2 2970.0 -o ./vtka_CH4_1024k_REF.txt "$InputDir/BlockInterface_1024k_CH4_REF_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 1790 -w 10.0 --bbox0 1790.0 --bbox1 1790.0 --bbox2 1790.0 -o ./vtka_C6H12_1024k.txt "$InputDir/BlockInterface_1024k_C6H12_node0_10000.vtu" )
#(cd $VTKA && ./vtk-analysis -n 1790 -w 10.0 --bbox0 1790.0 --bbox1 1790.0 --bbox2 1790.0 -o ./vtka_C6H12_1024k_REF.txt "$InputDir/BlockInterface_1024k_C6H12_REF_node0_10000.vtu" )
#echo "[Progress] 6/6"