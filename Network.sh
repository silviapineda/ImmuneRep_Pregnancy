#!/bin/bash
FILES=/Users/Pinedasans/ImmuneRep_Pregnancy/Results/Network/IGHM_unmutated/edges*
export PATH=/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin:$PATH
for f in $FILES
do  
    var=$(basename $f)
    echo $var
    echo $f
	java -cp nucleotides-assembly-1.0_second.jar  Graph $f $var.outcome.txt
done

