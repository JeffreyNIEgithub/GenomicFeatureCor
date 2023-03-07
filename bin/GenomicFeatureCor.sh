##!/bin/bash

##################################################################################################################
# 两两bed之间距离计算（bedtools）
# 20201210 
##################################################################################################################
#参数传递
while getopts "t:a:l:h" opt; do
  case $opt in
    t)
      targetlist=$OPTARG
      ;;
    a)
      arrowslist=$OPTARG
      ;;
    l)
      chrlength=$OPTARG
      ;;
    h)
     echo ""
      ;;
#    \?)
#     echo "Invalid option: -$OPTARG"
#     ;;
  esac
done

######################################################
#帮助文档
display_usage() {
        echo -e "\n\tThe script is designed to calculate the location information(bed file) \n"
                echo -e "\tUsage:  bash bed_dist.sh -t [] -a []
                echo -e "\tExample:bash bed_dist.sh -t targetlist -a arrowslist
        echo -e "\t-t: Necessary parameter. the target bed name\n"
        echo -e "\t-a: Necessary parameter. the bed name\n"
	echo -e "\t-l: Necessary parameter. the length of each chromosome\n"
                echo -e "\t-h or --help: Usage\n"
        }

# check whether user had supplied -h or --help . If yes display usage
        if [[ ( $* == "--help") ||  $* == "-h" ]]
        then
                display_usage
                exit 0
        fi

######################################################
#变量判定

if [ ! $targetlist ] || [ ! $arrowslist ] || [ ! $chrlength ];then
        echo -e "\n\tERROR: Missing necessary files."
        exit 1
fi

##############################################################################################################
#加载程序环境变量
pwd=`pwd`

bedtools=/home/.opt/biosoft/miniconda3_for_pb-assembly/bin/bedtools

####################################################################################################################
# 获得随机bed文件

for b in `cat ${arrowslist}`
do 
$bedtools shuffle -i $pwd/data/${b}.bed.txt  -g ${chrlength} -seed 100 |sort -k1,1 -k2,2n > $pwd/data/shuffle.${b}.bed
done

for a in `cat ${targetlist}`
do 
$bedtools shuffle -i $pwd/data/${a}.bed.txt  -g ${chrlength} -seed 100 |sort -k1,1 -k2,2n > $pwd/data/shuffle.${a}.bed
done

##################################################################################################################
#4. 距离数据准备

cd $pwd

##################################################################
#4.1. 绝对距离的计算
for i in `cat ${targetlist}`
do
mkdir ${i} && cd ${i}
for j in `cat $pwd/${arrowslist}`
do
##################################
#All
if [[ ${i} != ${j} ]];then
mkdir ${i}.vs.${j} && cd ${i}.vs.${j}
 $bedtools closest -a $pwd/data/${i}.bed.txt -b $pwd/data/${j}.bed.txt |awk '{if($5 != ".")print $0}'|$bedtools overlap -i stdin -cols 2,3,6,7 |awk -va="${i}.vs.${j}" -vb="All" -vc="Observed" '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,a,b,c,$9*-1}' > ${i}.vs.${j}.abs.All.txt
 $bedtools closest -a $pwd/data/${i}.bed.txt -b $pwd/data/shuffle.${j}.bed |awk '{if($5 != ".")print $0}'|$bedtools overlap -i stdin -cols 2,3,6,7 |awk -va="${i}.vs.Shuffle.${j}" -vb="All" -vc="Random" '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,a,b,c,$9*-1}' > ${i}.vs.${j}.shuffle.All.txt
 $bedtools closest -a $pwd/data/shuffle.${i}.bed -b $pwd/data/${j}.bed.txt|awk '{if($5 != ".")print $0}' |$bedtools overlap -i stdin -cols 2,3,6,7 |awk -va="${i}.Shuffle.vs.${j}" -vb="All" -vc="Random" '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,a,b,c,$9*-1}' > ${i}.shuffle.vs.${j}.All.txt
 
##################################
#upstream
 $bedtools closest -a $pwd/data/${i}.bed.txt -b $pwd/data/${j}.bed.txt  -D a |awk '{if($5 != ".")print $0}' |$bedtools sort|awk -va="${i}.vs.${j}" -vb="upstream" -vc="Observed" '{OFS="\t"}{if($9<0)print $1,$2,$3,$4,$5,$6,$7,$8,a,b,c,$9*-1}' > ${i}.vs.${j}.abs.upstream.txt
 $bedtools closest -a $pwd/data/${i}.bed.txt -b $pwd/data/shuffle.${j}.bed -D a |awk '{if($5 != ".")print $0}' |$bedtools sort|awk -va="${i}.vs.Shuffle.${j}" -vb="upstream" -vc="Random" '{OFS="\t"}{if($9<0)print $1,$2,$3,$4,$5,$6,$7,$8,a,b,c,$9*-1}' > ${i}.vs.${j}.shuffle.upstream.txt
 $bedtools closest -a $pwd/data/shuffle.${i}.bed  -b $pwd/data/${j}.bed.txt -D a |awk '{if($5 != ".")print $0}' |$bedtools sort|awk -va="${i}.Shuffle.vs.${j}" -vb="upstream" -vc="Random" '{OFS="\t"}{if($9<0)print $1,$2,$3,$4,$5,$6,$7,$8,a,b,c,$9*-1}' > ${i}.shuffle.vs.${j}.upstream.txt

##################################
#downstream
 $bedtools closest -a $pwd/data/${i}.bed.txt -b $pwd/data/${j}.bed.txt  -D a |awk '{if($5 != ".")print $0}' |awk -va="${i}.vs.${j}" -vb="downstream" -vc="Observed" '{OFS="\t"}{if ($9>0) print $1,$2,$3,$4,$5,$6,$7,$8,a,b,c,$9*-1}' > ${i}.vs.${j}.abs.downstream.txt
 $bedtools closest -a $pwd/data/${i}.bed.txt -b $pwd/data/shuffle.${j}.bed -D a |awk '{if($5 != ".")print $0}' |awk -va="${i}.vs.Shuffle.${j}" -vb="downstream" -vc="Random" '{OFS="\t"}{if($9>0)print $1,$2,$3,$4,$5,$6,$7,$8,a,b,c,$9*-1}' > ${i}.vs.${j}.shuffle.downstream.txt
 $bedtools closest -a $pwd/data/shuffle.${i}.bed -b $pwd/data/${j}.bed.txt -D a |awk '{if($5 != ".")print $0}' |awk -va="${i}.Shuffle.vs.${j}" -vb="downstream" -vc="Random" '{OFS="\t"}{if($9>0)print $1,$2,$3,$4,$5,$6,$7,$8,a,b,c,$9*-1}' > ${i}.shuffle.vs.${j}.downstream.txt

##################################
#overlap
 cat ${i}.vs.${j}.abs.All.txt |awk -va="${i}.vs.${j}" -vb="overlap" -vc="Observed" '{OFS="\t"}{if($12<=0)print $1,$2,$3,$4,$5,$6,$7,$8,a,b,c,$12*-1}' > ${i}.vs.${j}.overlap.txt
 cat ${i}.vs.${j}.shuffle.All.txt|awk -va="${i}.vs.Shuffle.${j}" -vb="overlap" -vc="Random" '{OFS="\t"}{if($12<=0)print $1,$2,$3,$4,$5,$6,$7,$8,a,b,c,$12*-1}' > ${i}.vs.${j}.shuffle.overlap.txt
 cat ${i}.shuffle.vs.${j}.All.txt|awk -va="${i}.Shuffle.vs.${j}" -vb="overlap" -vc="Random" '{OFS="\t"}{if($12<=0)print $1,$2,$3,$4,$5,$6,$7,$8,a,b,c,$12*-1}' > ${i}.shuffle.vs.${j}.overlap.txt

 ##################################
#合并

cat ${i}.vs.${j}.abs.All.txt ${i}.vs.${j}.shuffle.All.txt  ${i}.vs.${j}.abs.upstream.txt ${i}.vs.${j}.shuffle.upstream.txt  ${i}.vs.${j}.abs.downstream.txt ${i}.vs.${j}.shuffle.downstream.txt  ${i}.vs.${j}.overlap.txt ${i}.vs.${j}.shuffle.overlap.txt  > ${i}.vs.${j}.absdis.tsv

#取log的数据
awk '{if($12 <= 0)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"($12*-1)+1;else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12+1}' ${i}.vs.${j}.absdis.tsv  > ${i}.vs.${j}.absdis_log.tsv

#生成不取log的数据
cat ${i}.vs.${j}.absdis.tsv |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' > ${i}.vs.${j}.absdis_raw.tsv

##################################################################
#4.2.相对距离的计算
$bedtools reldist -a $pwd/data/${i}.bed.txt -b $pwd/data/${j}.bed.txt |awk -va="${i}.vs.${j}" -vb="Observed" '{OFS="\t"}NR>1{print $1,$2,$3,a,b}'> ${i}.vs.${j}.reldist.All.txt
$bedtools reldist -a $pwd/data/${i}.bed.txt -b $pwd/data/shuffle.${j}.bed |awk -va="${i}.vs.Shuffle.${j}" -vb="Random" '{OFS="\t"}NR>1{print $1,$2,$3,a,b}' > ${i}.vs.${j}.shuffle.reldist.All.txt
$bedtools reldist -a $pwd/data/shuffle.${i}.bed -b $pwd/data/${j}.bed.txt |awk -va="${i}.Shuffle.vs.${j}" -vb="Random" '{OFS="\t"}NR>1{print $1,$2,$3,a,b}' > ${i}.shuffle.vs.${j}.reldist.All.txt

##################################################################
#5.综合图
#数据准备-相对距离
cat ${i}.vs.${j}.reldist.All.txt ${i}.vs.${j}.shuffle.reldist.All.txt > ${i}.vs.${j}.reldist.All.tsv
cat ${i}.vs.${j}.reldist.All.txt ${i}.shuffle.vs.${j}.reldist.All.txt > ${i}.vs.${j}.reldist_v2.All.tsv
#数据准备-绝对距离
cat ${i}.vs.${j}.abs.All.txt ${i}.vs.${j}.abs.upstream.txt ${i}.vs.${j}.abs.downstream.txt ${i}.vs.${j}.overlap.txt ${i}.shuffle.vs.${j}.All.txt ${i}.shuffle.vs.${j}.upstream.txt ${i}.shuffle.vs.${j}.downstream.txt ${i}.shuffle.vs.${j}.overlap.txt > ${i}.vs.${j}.absdis_O_R_self.tsv
##overlap箱线图——all
cat ${i}.vs.${j}.overlap.txt ${i}.vs.${j}.shuffle.overlap.txt > overlap_boxplot.tsv


rm *txt

Rscript SinglePlot.r ./ ${i}.vs.${j}
cd ..
fi
done
cd $pwd/${i}
cat ./*/*.reldist.All.tsv |sed 's/Shuffle.//g' > ${i}_others.reldist.All.tsv
cat ./*/*.absdis_O_R_self.tsv |grep 'overlap' > ${i}_others_overlap_absdis_O_R_self.tsv

awk '{print $1"_"$2"_"$3"_"$4"\t"$9"\t"$10"\t"$11}' ${i}_others_overlap_absdis_O_R_self.tsv |grep 'Observed'|sort|uniq|awk '{print $2}' |sort|uniq -c |awk -va="Observed" '{OFS="\t"}{print $2"\t"$1"\t"a}' > overlap_num_Observed
awk '{print $1"_"$2"_"$3"_"$4"\t"$9"\t"$10"\t"$11}' ${i}_others_overlap_absdis_O_R_self.tsv |grep 'Random'|sort|uniq|awk '{print $2}' |sort|uniq -c |awk -va="Random" '{OFS="\t"}{print $2"\t"$1"\t"a}' > overlap_num_Random
cat overlap_num_Observed overlap_num_Random |sed 's/Shuffle.//g'|sed 's/${i}.vs.//g'|sed '1iType\tnumber\tGroup' > overlap_number
cat ./*/overlap_boxplot.tsv |sed 's/Shuffle.//g' > overlap_boxplot_v2.tsv
Rscript CombinedPlot.r ./ ${i}.vs.${j} ${i}
rm overlap_num_Observed overlap_num_Random
cd $pwd
done




