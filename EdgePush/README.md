# VLDB 2022: Edge-based Local Push for Personalized PageRank



## Tested Environment:
- Ubuntu 16.04.10
- g++ 11
- GCC 5.4.0


## Datasets:
We provide the data sources in the paper. Notably, after downloading the raw datasets, we first need to convert them to motif-based weighted graphs. In all of our experiments, we employ the codes released by MAPPR [Yin et al, KDD2017] to convert the raw datasets to motif-based weighted graphs. The source code of MAPPR is available at "http://snap.stanford.edu/mappr/". Please save the obtained weighted graphs in "EdgePush-master/dataset/" and rename them as "youtube.txt, soc-LiveJournal1.txt, indochina-2004.txt or orkut-links.txt", separately. 


[Yin et al, KDD2017] Yin H, Benson A R, Leskovec J, et al. Local higher-order graph clustering[C]//Proceedings of the 23rd ACM SIGKDD international conference on knowledge discovery and data mining. 2017: 555-564.




## Commands:
We include the fundamental commands in the script file: "EdgePush-master/run_script.sh". To execute our codes automatically, please use the following bash commands: 
```
bash run_script.sh
```


Alternatively, our codes can be executed mannually. Specifically, to compile the codes: 
```
rm SSPPR
make
```

To run powermethod for groundtruths: 
```
./SSPPR -d ./dataset/grqc/ -f grqc -algo powermethod -qn 10 -a 0.2
```


To run EdgePush: 
```
./SSPPR -d ./dataset/grqc/ -f grqc -algo EdgePush -qn 10 -a 0.2
```

Notably, in our experiments, we use the public codes of MAPPR, FORA, SpeedPPR, PowForPush as our baselines, which are available from: 
* MAPPR: http://snap.stanford.edu/mappr/
* FORA: https://github.com/wangsibovictor/fora
* SpeedPPR, PowForPush: https://github.com/wuhao-wu-jiang/Personalized-PageRank




## Parameters:
- -d \<directory\> 
- -f \<filelabel\>
- -algo \<algorithm\>
- [-e \<rmax\> (0.1 by default)]
- [-qn \<the number of query nodes\> (10 by default)]
- [-a \<alpha\> (0.2 by default)]



## Remark:
* EdgePush-master/datatset/: containing the datasets 
* EdgePush-master/query/: containing the files for query nodes
* EdgePush-master/result/: containing the SSPPR approximation. 
