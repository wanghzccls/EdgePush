para1_arr=("0.001" "0.0001" "1e-05" "1e-06" "1e-07" "1e-08")
para2_arr=(0.2)


data_arr=("grqc") # a small dataset
#data_arr=("youtube")
#data_arr=("youtube" "soc-LiveJournal1" "indochina-2004" "orkut-links") # four real graphs used in our paper

echo "=== SSPPR Approximation ==="
rm SSPPR
make

for data_name in "${data_arr[@]}"
do	
	for((j=0;j<${#para2_arr[@]};j++))
	do
		echo "./SSPPR -d ./dataset/${data_name}/ -f ${data_name} -algo powermethod -qn 10 -a ${para2_arr[$j]}" |bash;
		for((i=0;i<${#para1_arr[@]};i++))
		do
			echo "./SSPPR -d ./dataset/${data_name}/ -f ${data_name} -algo EdgePush -e ${para1_arr[$i]} -qn 10 -a ${para2_arr[$j]}" |bash;
		done
	done
done




