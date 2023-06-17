#!/usr/bin/zsh

MarDynDir=../../../cmake-build-debug/src/
InputDir=$(pwd)

## 1k simulations
#echo "Running 1k simulations"
#echo "[Progress] 0/24"
#(cd $MarDynDir && rm -r -f BI1k_C6H12.log)
#rm -r -f tmp.conf.xml
#sed -r 's#<output>#<!--output>#g' "$InputDir/1k/config_C6H12.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
#(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI1k_C6H12.log)
#echo "[Progress] 1/24"
#(cd $MarDynDir && rm -r -f BI1k_C6H12_REF.log)
#rm -r -f tmp.conf.xml
#sed -r 's#<output>#<!--output>#g' "$InputDir/1k/config_C6H12_REF.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
#(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI1k_C6H12_REF.log)
#
#echo "[Progress] 2/24"
#(cd $MarDynDir && rm -r -f BI1k_CH4.log)
#rm -r -f tmp.conf.xml
#sed -r 's#<output>#<!--output>#g' "$InputDir/1k/config_CH4.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
#(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI1k_CH4.log)
#echo "[Progress] 3/24"
#(cd $MarDynDir && rm -r -f BI1k_CH4_REF.log)
#rm -r -f tmp.conf.xml
#sed -r 's#<output>#<!--output>#g' "$InputDir/1k/config_CH4_REF.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
#(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI1k_CH4_REF.log)
#
## 4k simulations
#echo "Running 4k simulations"
#echo "[Progress] 4/24"
#(cd $MarDynDir && rm -r -f BI4k_C6H12.log)
#rm -r -f tmp.conf.xml
#sed -r 's#<output>#<!--output>#g' "$InputDir/4k/config_C6H12.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
#(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI4k_C6H12.log)
#echo "[Progress] 5/24"
#(cd $MarDynDir && rm -r -f BI4k_C6H12_REF.log)
#rm -r -f tmp.conf.xml
#sed -r 's#<output>#<!--output>#g' "$InputDir/4k/config_C6H12_REF.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
#(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI4k_C6H12_REF.log)
#
#echo "[Progress] 6/24"
#(cd $MarDynDir && rm -r -f BI4k_CH4.log)
#rm -r -f tmp.conf.xml
#sed -r 's#<output>#<!--output>#g' "$InputDir/4k/config_CH4.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
#(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI4k_CH4.log)
#echo "[Progress] 7/24"
#(cd $MarDynDir && rm -r -f BI4k_CH4_REF.log)
#rm -r -f tmp.conf.xml
#sed -r 's#<output>#<!--output>#g' "$InputDir/4k/config_CH4_REF.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
#(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI4k_CH4_REF.log)
#
## 16k simulations
#echo "Running 16k simulations"
#echo "[Progress] 8/24"
#(cd $MarDynDir && rm -r -f BI16k_C6H12.log)
#rm -r -f tmp.conf.xml
#sed -r 's#<output>#<!--output>#g' "$InputDir/16k/config_C6H12.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
#(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI16k_C6H12.log)
#echo "[Progress] 9/24"
#(cd $MarDynDir && rm -r -f BI16k_C6H12_REF.log)
#rm -r -f tmp.conf.xml
#sed -r 's#<output>#<!--output>#g' "$InputDir/16k/config_C6H12_REF.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
#(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI16k_C6H12_REF.log)
#
#echo "[Progress] 10/24"
#(cd $MarDynDir && rm -r -f BI16k_CH4.log)
#rm -r -f tmp.conf.xml
#sed -r 's#<output>#<!--output>#g' "$InputDir/16k/config_CH4.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
#(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI16k_CH4.log)
#echo "[Progress] 11/24"
#(cd $MarDynDir && rm -r -f BI16k_CH4_REF.log)
#rm -r -f tmp.conf.xml
#sed -r 's#<output>#<!--output>#g' "$InputDir/16k/config_CH4_REF.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
#(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI16k_CH4_REF.log)
#
## 64k simulations
#echo "Running 64k simulations"
#echo "[Progress] 12/24"
#(cd $MarDynDir && rm -r -f BI64k_C6H12.log)
#rm -r -f tmp.conf.xml
#sed -r 's#<output>#<!--output>#g' "$InputDir/64k/config_C6H12.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
#(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI64k_C6H12.log)
#echo "[Progress] 13/24"
#(cd $MarDynDir && rm -r -f BI64k_C6H12_REF.log)
#rm -r -f tmp.conf.xml
#sed -r 's#<output>#<!--output>#g' "$InputDir/64k/config_C6H12_REF.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
#(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI64k_C6H12_REF.log)

echo "[Progress] 14/24"
(cd $MarDynDir && rm -r -f BI64k_CH4.log)
rm -r -f tmp.conf.xml
sed -r 's#<output>#<!--output>#g' "$InputDir/64k/config_CH4.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI64k_CH4.log)
echo "[Progress] 15/24"
(cd $MarDynDir && rm -r -f BI64k_CH4_REF.log)
rm -r -f tmp.conf.xml
sed -r 's#<output>#<!--output>#g' "$InputDir/64k/config_CH4_REF.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI64k_CH4_REF.log)

# 256k simulations
echo "Running 256k simulations"
echo "[Progress] 16/24"
(cd $MarDynDir && rm -r -f BI256k_C6H12.log)
rm -r -f tmp.conf.xml
sed -r 's#<output>#<!--output>#g' "$InputDir/256k/config_C6H12.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI256k_C6H12.log)
echo "[Progress] 17/24"
(cd $MarDynDir && rm -r -f BI256k_C6H12_REF.log)
rm -r -f tmp.conf.xml
sed -r 's#<output>#<!--output>#g' "$InputDir/256k/config_C6H12_REF.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI256k_C6H12_REF.log)

echo "[Progress] 18/24"
(cd $MarDynDir && rm -r -f BI256k_CH4.log)
rm -r -f tmp.conf.xml
sed -r 's#<output>#<!--output>#g' "$InputDir/256k/config_CH4.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI256k_CH4.log)
echo "[Progress] 19/24"
(cd $MarDynDir && rm -r -f BI256k_CH4_REF.log)
rm -r -f tmp.conf.xml
sed -r 's#<output>#<!--output>#g' "$InputDir/256k/config_CH4_REF.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI256k_CH4_REF.log)

# 1024k simulations
echo "Running 1024k simulations"
echo "[Progress] 20/24"
(cd $MarDynDir && rm -r -f BI1024k_C6H12.log)
rm -r -f tmp.conf.xml
sed -r 's#<output>#<!--output>#g' "$InputDir/1024k/config_C6H12.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI1024k_C6H12.log)
echo "[Progress] 21/24"
(cd $MarDynDir && rm -r -f BI1024k_C6H12_REF.log)
rm -r -f tmp.conf.xml
sed -r 's#<output>#<!--output>#g' "$InputDir/1024k/config_C6H12_REF.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI1024k_C6H12_REF.log)

echo "[Progress] 22/24"
(cd $MarDynDir && rm -r -f BI1024k_CH4.log)
rm -r -f tmp.conf.xml
sed -r 's#<output>#<!--output>#g' "$InputDir/1024k/config_CH4.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI1024k_CH4.log)
echo "[Progress] 23/24"
(cd $MarDynDir && rm -r -f BI1024k_CH4_REF.log)
rm -r -f tmp.conf.xml
sed -r 's#<output>#<!--output>#g' "$InputDir/1024k/config_CH4_REF.xml" | sed -r 's#</output>#<output-->#g' - | sed -r 's#\./\.\./#\./#g' -  > tmp.conf.xml
(cd $MarDynDir && ./MarDyn "$InputDir/tmp.conf.xml" >> BI1024k_CH4_REF.log)
echo "[Progress] 24/24"