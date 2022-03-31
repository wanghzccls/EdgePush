# VLDB 2022: Edge-based Local Push for Personalized PageRank

## Citation:
```
@article{wang2022edge,
  title={Edge-based Local Push for Personalized PageRank},
  author={Wang, Hanzhi and Wei, Zhewei and Gan, Junhao and Yuan, Ye and Du, Xiaoyong and Wen, Ji-Rong},
  journal={arXiv preprint arXiv:2203.07937},
  year={2022}
}
```



## Tested Environment:
- Ubuntu 16.04.10
- g++ 11
- GCC 5.4.0


## Datasets:
To conduct the experiments on motif-based weighted graphs, we need to convert the original unweighted graphs to motif-based weighted graphs using the codes released by MAPPR [Yin et al, KDD2017] ("http://snap.stanford.edu/mappr/"). Please save the obtained weighted graphs in "./dataset/" and rename them as "youtube.txt, soc-LiveJournal1.txt, indochina-2004.txt or orkut-links.txt", separately. 


[Yin et al, KDD2017] Yin H, Benson A R, Leskovec J, et al. Local higher-order graph clustering[C]//Proceedings of the 23rd ACM SIGKDD international conference on knowledge discovery and data mining. 2017: 555-564.




## Commands:
Please execute the script ("./run_script.sh") to run our codes automatically: 
```
bash run_script.sh
```


Alternatively, the codes can be compiled and executed mannually as following: 

To compile the codes: 
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

Notably, in our experiments, we employ MAPPR, FORA, SpeedPPR, PowForPush as our baselines, which are available at: 
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
* "./datatset/" for the datasets 
* "./query/" for the query nodes
* "./result/" for the SSPPR approximation results. 
