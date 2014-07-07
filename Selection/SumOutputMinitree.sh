#!/bin/bash
#list=`find -name "v05_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_START42_V11_closure*.out"`
syntax="Syntax ${0} {sampleName}"
if [[ -z ${1} ]]
then
  echo ${syntax}
#  exit 1
fi

list=`find -name "${1}*.out"`

a=a2=b=b2=c=c2=d=d2=e=e2=l3=l4=m=m2=n=n2=o=o2=p=p2=q=q2=r1=r2=s=s2=t=t2=u=u2=v=v2=w=w2=x=x2=y=y2=z=z2=alpha=alpha2=beta=beta2=gamma=gamma2=delta=delta2=0;

for fichier in $list 
do
	echo "$fichier"
	let "a+=`grep "TOTALnbMuonsAfterID\[0\]" $fichier | cut -f 2`"
	let "a2+=`grep "TOTALnbMuonsAfterID\[0\]" $fichier | cut -f 5`"
	let "b+=`grep "TOTALnbMuonsAfterID\[1\]" $fichier | cut -f 2`"
        let "b2+=`grep "TOTALnbMuonsAfterID\[1\]" $fichier | cut -f 5`"
	let "c+=`grep "TOTALnbMuonsAfterID\[2\]" $fichier | cut -f 2`"
        let "c2+=`grep "TOTALnbMuonsAfterID\[2\]" $fichier | cut -f 5`"
	let "d+=`grep "TOTALnbMuonsAfterID\[3\]" $fichier | cut -f 2`"
        let "d2+=`grep "TOTALnbMuonsAfterID\[3\]" $fichier | cut -f 5`"
	let "e+=`grep "TOTALnbMuonsAfterID\[4\]" $fichier | cut -f 2`"
        let "e2+=`grep "TOTALnbMuonsAfterID\[4\]" $fichier | cut -f 5`"

	let "l3+=`grep "TOTALnbDimuonsAfterID\[0\]" $fichier | cut -f 2`"
        let "l4+=`grep "TOTALnbEventsAfterDimuonID\[0\]" $fichier | cut -f 5`"
	let "m+=`grep "TOTALnbDimuonsAfterID\[1\]" $fichier | cut -f 2`"
        let "m2+=`grep "TOTALnbEventsAfterDimuonID\[1\]" $fichier | cut -f 5`"
	let "n+=`grep "TOTALnbDimuonsAfterID\[2\]" $fichier | cut -f 2`"
        let "n2+=`grep "TOTALnbEventsAfterDimuonID\[2\]" $fichier | cut -f 5`"

	let "o+=`grep "TOTALnbPhotonsAfterID\[0\]" $fichier | cut -f 2`"
        let "o2+=`grep "TOTALnbEventsAfterPhotonID\[0\]" $fichier | cut -f 5`"
	let "p+=`grep "TOTALnbPhotonsAfterID\[1\]" $fichier | cut -f 2`"
        let "p2+=`grep "TOTALnbEventsAfterPhotonID\[1\]" $fichier | cut -f 5`"
	let "q+=`grep "TOTALnbPhotonsAfterID\[2\]" $fichier | cut -f 2`"
        let "q2+=`grep "TOTALnbEventsAfterPhotonID\[2\]" $fichier | cut -f 5`"
	let "r1+=`grep "TOTALnbPhotonsAfterID\[3\]" $fichier | cut -f 2`"
        let "r2+=`grep "TOTALnbEventsAfterPhotonID\[3\]" $fichier | cut -f 5`"
	let "s+=`grep "TOTALnbPhotonsAfterID\[4\]" $fichier | cut -f 2`"
        let "s2+=`grep "TOTALnbEventsAfterPhotonID\[4\]" $fichier | cut -f 5`"
	let "t+=`grep "TOTALnbPhotonsAfterID\[5\]" $fichier | cut -f 2`"
        let "t2+=`grep "TOTALnbEventsAfterPhotonID\[5\]" $fichier | cut -f 5`"

	let "u+=`grep "TOTALnbMuMuGammaAfterID\[0\]" $fichier | awk '{print $3}'`"
        let "u2+=`grep "TOTALnbEventsAfterMuMuGammaID\[0\]" $fichier | awk '{print $6}'`"
	let "v+=`grep "TOTALnbMuMuGammaAfterID\[1\]" $fichier | awk '{print $3}'`"
        let "v2+=`grep "TOTALnbEventsAfterMuMuGammaID\[1\]" $fichier | awk '{print $6}'`"
	let "w+=`grep "TOTALnbMuMuGammaAfterID\[2\]" $fichier | awk '{print $3}'`"
        let "w2+=`grep "TOTALnbEventsAfterMuMuGammaID\[2\]" $fichier | awk '{print $6}'`"
	let "x+=`grep "TOTALnbMuMuGammaAfterID\[3\]" $fichier | awk '{print $3}'`"
        let "x2+=`grep "TOTALnbEventsAfterMuMuGammaID\[3\]" $fichier | awk '{print $6}'`"
	let "y+=`grep "TOTALnbMuMuGammaAfterID\[4\]" $fichier | awk '{print $3}'`"
        let "y2+=`grep "TOTALnbEventsAfterMuMuGammaID\[4\]" $fichier | awk '{print $6}'`"
	let "z+=`grep "TOTALnbMuMuGammaAfterID\[5\]" $fichier | awk '{print $3}'`"
        let "z2+=`grep "TOTALnbEventsAfterMuMuGammaID\[5\]" $fichier | awk '{print $6}'`"
	let "alpha+=`grep "TOTALnbMuMuGammaAfterID\[6\]" $fichier | awk '{print $3}'`"
        let "alpha2+=`grep "TOTALnbEventsAfterMuMuGammaID\[6\]" $fichier | awk '{print $6}'`"

	let "beta+=`grep "TOTALnbMuMuGammaAfterID\[7\]" $fichier | awk '{print $3}'`"
        let "beta2+=`grep "TOTALnbEventsAfterMuMuGammaID\[7\]" $fichier | awk '{print $6}'`"

	let "gamma+=`grep "TOTALnbMuMuGammaAfterID\[8\]" $fichier | awk '{print $3}'`"
        let "gamma2+=`grep "TOTALnbEventsAfterMuMuGammaID\[8\]" $fichier | awk '{print $6}'`"

	let "delta+=`grep "TOTALnbMuMuGammaAfterID\[9\]" $fichier | awk '{print $3}'`"
        let "delta2+=`grep "TOTALnbEventsAfterMuMuGammaID\[9\]" $fichier | awk '{print $6}'`"



done  

bper=$(echo "scale=3;$b/$a * 100" | bc -l)
b2per=$(echo "scale=3;$b2/$a2 * 100" | bc -l)

cper=$(echo "scale=3;$c/$b * 100" | bc -l)
c2per=$(echo "scale=3;$c2/$b2 * 100" | bc -l)

dper=$(echo "scale=3;$d/$c * 100" | bc -l)
d2per=$(echo "scale=3;$e2/$c2 * 100" | bc -l)

eper=$(echo "scale=3;$e/$d * 100" | bc -l)
e2per=$(echo "scale=3;$e2/$d2 * 100" | bc -l)


l3per=$(echo "scale=3;$l3/$e * 100" | bc -l)
l4per=$(echo "scale=3;$l4/$e2 * 100" | bc -l)

mper=$(echo "scale=3;$m/$l3 * 100" | bc -l)
m2per=$(echo "scale=3;$m2/$l4 * 100" | bc -l)

nper=$(echo "scale=3;$n/$m * 100" | bc -l)
n2per=$(echo "scale=3;$n2/$m2 * 100" | bc -l)



oper=$(echo "scale=3;$o/$n * 100" | bc -l)
o2per=$(echo "scale=3;$o2/$n2 * 100" | bc -l)

pper=$(echo "scale=3;$p/$o * 100" | bc -l)
p2per=$(echo "scale=3;$p2/$o2 * 100" | bc -l)

qper=$(echo "scale=3;$q/$p * 100" | bc -l)
q2per=$(echo "scale=3;$q2/$p2 * 100" | bc -l)

r1per=$(echo "scale=3;$r1/$q * 100" | bc -l)
r2per=$(echo "scale=3;$r2/$q2 * 100" | bc -l)

sper=$(echo "scale=3;$s/$r1 * 100" | bc -l)
s2per=$(echo "scale=3;$s2/$r2 * 100" | bc -l)

tper=$(echo "scale=3;$t/$s * 100" | bc -l)
t2per=$(echo "scale=3;$t2/$s2 * 100" | bc -l)




uper=$(echo "scale=3;$u/$t * 100" | bc -l)
u2per=$(echo "scale=3;$u2/$t2 * 100" | bc -l)

vper=$(echo "scale=3;$v/$u * 100" | bc -l)
v2per=$(echo "scale=3;$v2/$u2 * 100" | bc -l)

wper=$(echo "scale=3;$w/$v * 100" | bc -l)
w2per=$(echo "scale=3;$w2/$v2 * 100" | bc -l)

xper=$(echo "scale=3;$x/$w * 100" | bc -l)
x2per=$(echo "scale=3;$x2/$w2 * 100" | bc -l)

yper=$(echo "scale=3;$y/$x * 100" | bc -l)
y2per=$(echo "scale=3;$y2/$x2 * 100" | bc -l)

zper=$(echo "scale=3;$z/$y * 100" | bc -l)
z2per=$(echo "scale=3;$z2/$y2 * 100" | bc -l)

alphaper=$(echo "scale=3;$alpha/$z * 100" | bc -l)
alpha2per=$(echo "scale=3;$alpha2/$z2 * 100" | bc -l)

betaper=$(echo "scale=3;$beta/$alpha * 100" | bc -l)
beta2per=$(echo "scale=3;$beta2/$alpha2 * 100" | bc -l)

gammaper=$(echo "scale=3;$gamma/$beta * 100" | bc -l)
gamma2per=$(echo "scale=3;$gamma2/$beta2 * 100" | bc -l)

deltaper=$(echo "scale=3;$delta/$gamma * 100" | bc -l)
delta2per=$(echo "scale=3;$delta2/$gamma2 * 100" | bc -l)

echo "TOTALnbMuonsAfterID[0]= ${a}    TOTALnbEventsAfterMuonID[0]= ${a2}"
echo "TOTALnbMuonsAfterID[1]= ${b}, ${bper}   TOTALnbEventsAfterMuonID[1]= ${b2}, ${b2per}"
echo "TOTALnbMuonsAfterID[2]= ${c}, ${cper}    TOTALnbEventsAfterMuonID[2]= ${c2}, ${c2per}"
echo "TOTALnbMuonsAfterID[3]= ${d}, ${dper}    TOTALnbEventsAfterMuonID[3]= ${d2}, ${d2per}"
echo "TOTALnbMuonsAfterID[4]= ${e}, ${eper}    TOTALnbEventsAfterMuonID[4]= ${e2}, ${e2per}"
echo ""
echo "TOTALnbDimuonsAfterID[0]= ${l3},     TOTALnbEventsAfterDimuonID[0]= ${l4}, ${l4per}"
echo "TOTALnbDimuonsAfterID[1]= ${m}, ${mper}    TOTALnbEventsAfterDimuonID[1]= ${m2}, ${m2per}"
echo "TOTALnbDimuonsAfterID[2]= ${n}, ${nper}    TOTALnbEventsAfterDimuonID[2]= ${n2}, ${n2per}"
echo ""
echo "TOTALnbPhotonsAfterID[0]= ${o},    TOTALnbEventsAfterPhotonID[0]= ${o2}, ${o2per}"
echo "TOTALnbPhotonsAfterID[1]= ${p}, ${pper}    TOTALnbEventsAfterPhotonID[1]= ${p2}, ${p2per}"
echo "TOTALnbPhotonsAfterID[2]= ${q}, ${qper}    TOTALnbEventsAfterPhotonID[2]= ${q2}, ${q2per}"
echo "TOTALnbPhotonsAfterID[3]= ${r1}, ${r1per}    TOTALnbEventsAfterPhotonID[3]= ${r2}, ${r2per}"
echo "TOTALnbPhotonsAfterID[4]= ${s}, ${sper}    TOTALnbEventsAfterPhotonID[4]= ${s2}, ${s2per}"
echo "TOTALnbPhotonsAfterID[5]= ${t}, ${tper}    TOTALnbEventsAfterPhotonID[5]= ${t2}, ${t2per}"
echo ""
echo "TOTALnbMuMuGammaAfterID[0]= ${u},     TOTALnbEventsAfterMuMuGammaID[0]= ${u2}, ${u2per}"
echo "TOTALnbMuMuGammaAfterID[1]= ${v}, ${vper}    TOTALnbEventsAfterMuMuGammaID[1]= ${v2}, ${v2per}"
echo "TOTALnbMuMuGammaAfterID[2]= ${w}, ${wper}    TOTALnbEventsAfterMuMuGammaID[2]= ${w2}, ${w2per}"
echo "TOTALnbMuMuGammaAfterID[3]= ${x}, ${xper}    TOTALnbEventsAfterMuMuGammaID[3]= ${x2}, ${x2per}"
echo "TOTALnbMuMuGammaAfterID[4]= ${y}, ${yper}    TOTALnbEventsAfterMuMuGammaID[4]= ${y2}, ${y2per}"
echo "TOTALnbMuMuGammaAfterID[5]= ${z}, ${zper}    TOTALnbEventsAfterMuMuGammaID[5]= ${z2}, ${z2per}"
echo "TOTALnbMuMuGammaAfterID[6]= ${alpha}, ${alphaper}    TOTALnbEventsAfterMuMuGammaID[6]= ${alpha2}, ${alpha2per}"
echo ""
echo "TOTALnbMuMuGammaAfterID[7]= ${beta}, ${betaper}    TOTALnbEventsAfterMuMuGammaID[7]= ${beta2}, ${beta2per}"
echo "TOTALnbMuMuGammaAfterID[8]= ${gamma}, ${gammaper}    TOTALnbEventsAfterMuMuGammaID[8]= ${gamma2}, ${gamma2per}"
echo "TOTALnbMuMuGammaAfterID[9]= ${delta}, ${deltaper}    TOTALnbEventsAfterMuMuGammaID[9]= ${delta2}, ${delta2per}"



