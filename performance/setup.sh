#!/bin/bash
cubit="cubit -nojournal -nographics"
decomp="decomp --processors "
prefix=""

function mesh() {
i=$1
if [ $i == 0 ]; then
cp ${prefix}.jou.tmpl ${prefix}.jou
else
sed "s/{}/$i/g" ${prefix}.jou.tmpl > ${prefix}.jou
fi
$cubit ${prefix}.jou &> cubit.log
if [ ! -e ${prefix}.g ]; then
  echo "!!! cubit FAILED !!!"
  exit
fi
mv ${prefix}.g ${prefix}_$i.g
}

function split() {
i=$1
j=$2
$decomp $j ${prefix}_${i}.g &> decomp.log
if [ ! -e ${prefix}_${i}.g.${j}.0 ] && [ ! -e ${prefix}_${i}.g.${j}.00 ] ; then
  echo "!!! decomp FAILED !!!"
  exit
fi
}

function setup() {
prefix=$1
for i in 0 1 2 ; do
  echo "... creating mesh $i"
  mesh $i
  for j in 4 16 64 ; do
    split $i $j
  done
done
}

for d in `ls -d *_contact` ; do
  echo ">> setting up $d"
  cd $d
  setup $d
  cd ..
done
