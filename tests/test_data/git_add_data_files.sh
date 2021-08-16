# add relavent files for testing

ls -1 *.cif | gawk '$0!~/FINAL_/' | while read cif; do
    echo $cif

    nebrs=${cif}.dat.total.txt
    charge=Final_Charge_${cif}.txt
    zero=Final_Charge_Zeroed_${cif}.txt
    shift=Shift_Charge_${cif}.txt
    final=FINAL_${cif}
    git add $cif $nebr $charge $zero $shift $final

done
