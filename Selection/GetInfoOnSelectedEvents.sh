#!/bin/bash
# to get the event information on selected events & put it in a nice format
# Written by O. Bondu (September 2010)

syntax="$0 {selection}"
if [[ -z $1 ]]
then
	echo ${syntax}
	exit 1
fi

selection=$1

echo "*** LaTeXTables/Selected_${selection}.tex ***"
echo "*** Selected_${selection}.txt ***"
rm Selected_${selection}.txt

selectedOutfilesList=""
for outfiles in `'ls' Selected/${selection}/Run2010*/*out`
#for outfiles in `echo "Selected/${selection}/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/*.out"`
do
	nSelected=`grep "nSelected=" ${outfiles} | awk '{print $2}'`
	if [[ ${nSelected} -gt 0 ]]
	then
		selectedOutfilesList=`echo "${selectedOutfilesList}${outfiles} "`
	fi
done

echo "
\setlength{\tabcolsep}{3pt}
\begin{tabular}{ccc|cccccccccc}
\hline
\rowcolor[gray]{.9}run & lumi & event & $ m_{\mu\mu\gamma}$ & $ {E_T}_{\gamma}$ & $ {p_T}_{\mu_{close}}$ & $ {p_T}_{\mu_{far}}$ & $\Delta R_{(\gamma, \mu_{close})}$ & $\Delta R_{(\gamma, \mu_{far})}$ & $ m_{\mu\mu}$ & $ \eta_{\mu_{close}}$ & $ \eta_{\mu_{far}}$ & $ \eta_{\gamma}$
\\\\\rowcolor[gray]{.9} & & & $\giga\eV$ & $\giga\eV$ & $\giga\eV$ & $\giga\eV$ & & & $\giga\eV$ & & & 
\\\\\hline
" > LaTeXTables/Selected_${selection}.tex

for file in `echo "${selectedOutfilesList}"`
do
	numbers=`grep -n -w "\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*" ${file} | cut -d : -f -1`
	firstline=`echo "${numbers}" | head -n 1`
	lastline=`echo "${numbers}" | tail -n 1`
	let firstline=${firstline=}+3
	let lastline=${lastline=}-2
	sed -n "${firstline},${lastline} p" ${file} >> Selected_${selection}.txt
	sed -n "${firstline},${lastline} p" ${file} | sed -e 's/\t\t/$ \& \$/g' -e 's/$/\$\\\\/g' -e 's/^/\$/g' >> LaTeXTables/Selected_${selection}.tex
done

echo "\hline
\end{tabular}"  >> LaTeXTables/Selected_${selection}.tex


exit 0
