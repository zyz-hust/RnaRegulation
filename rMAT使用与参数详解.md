## 简介
rMAT是一款对RNA-seq数据进行差异可变剪切分析的软件。通过rMATS统计模型对不同样本进行可变剪切事件的表达定量，然后以likelihood-ratio  test计算P value 来表示两组样品在IncLevel水平上的差异，IncLevel利用Benjamini Hochberg算法对p value 校正得FDR值。rMATS可识别的可变剪切有五种：

- 1. SE ：skipped exon  外显子跳跃
- 2. A5SS ：alternative 5’ splice site 第一个外显子可变剪切
- 3. A3SS ：alternative 3’ splice site 最后一个外显子可变剪切
- 4. MXE ：mutually exclusive exons 外显子选择性跳跃
- 5. RI ：retained intron  内含子滞留

![Alternative Splicing Events.png](https://i.loli.net/2021/01/18/lhuzdS51YPeNtI3.png)

## 软件下载

```linux
# 仅官网或者下述网站进行下载
https://sourceforge.net/projects/rnaseq\-mats/files/MATS/rMATS.4.0.1.tgz/download

pip install numpy
sudo apt-get install libblas-dev liblas-dev liblapack-dev
sudo apt\-get install gfortran
```

或者使用conda 下载

```linux
conda install -c bioconda rmats
```

新版rMATS解包后会有两个rmat.py执行脚本，因此需要首先测试系统支持的那种
进入python 输入以下命令：

```python
import sys
print sys.maxunicode
# 如果出现1114111则说明需要使用rMATS-turbo-Linux-UCS4文件下的rmats.py;出现65535则说明需要使用rMATS-turbo-Linux-UCS2文件下rmats.py
```

## rMATS 参数设置
```linux
# 完整命令
/BioII/lulab\_b/chenyinghui/software/conda2/bin/python2 \
/BioII/lulab\_b/chenyinghui/software/rMATS/rmats\_turbo\_v4\_1\_0\_python\_2\_7/rmats.py \
--b1 /BioII/lulab\_b/liuxiaofan/project/ribo-shape/result\_liulab\_batch4/splicing\_lxf/wt.noUVB.txt \
--b2 /BioII/lulab\_b/liuxiaofan/project/ribo-shape/result\_liulab\_batch4/splicing\_lxf/wt.UVB.txt \
--gtf /BioII/lulab\_b/liuxiaofan/database/ATH/GTF/Arabidopsis\_thaliana.TAIR10.34.gtf \
--od /BioII/lulab\_b/liuxiaofan/project/ribo-shape/result\_liulab\_batch4/splicing\_lxf/wt.UVB-vs-noUVB \
-t paired \ 
--readLength 151 \ 
--cstat 0.0001  \
--tmp /BioII/lulab\_b/liuxiaofan/project/ribo-shape/result\_liulab\_batch4/splicing\_lxf/wt.UVB-vs-noUVB/tmp \
--nthread 4  \
--variable-read-length
```

* --b1 b1.txt 输入sample1的txt格式的文件。文件内以逗号分隔重复样本的bam文件名
* --b2 b2.txt 同--b1
* -gtf gtfFILE 需要输入的gtf文件
* --od outDir 所有输出文件的路径文件夹
* --nthread 设置线程数
* --cstat 设置splicing difference的阈值
* --readLength 151  给定读段长度
* --tmp  tmpDir 设置临时文件夹
* --variable-read-length 

## 结果文件
rMATS的结果文件是以各个可变剪切事件的分布。
XX 指代SE\RI\MXE\A5SS\A3SS 五项可变剪切时间。
* 1. fromGTF.XX.txt 系列： 直接从GTF文件和数据文件读出的结果
* 2. fromGTF.novelEvents.XX.txt : 从数据文件发现的新的可变剪切事件
* 3. XX.MATS.JC.txt 和XX.MATS.JECE.txt,是JC.raw.input.XX.txt 和 JCEC.raw.input.XX.txt 经过统计模型分析后的结果，多了P值和FDR值。
* 4. JC和JCEC的区别在于前者考虑跨越剪切位点的reads，而后者不仅考虑前者的reads还考虑到只比对到第一张图的条纹区(没有跨越剪切位点的reads)，一般只是用JC就够了。

**XX.MATS.JC.txt中包含的信息比较多，以SE.MATS.JC.txt为例：**
* 1. 1-5列分别为：ID、GeneID、geneSymbol、chr、strand
* 2. 6-11列分别为外显子位置信息：分别为exonStart_0base，exonEnd，upstreamES，upstreamEE，downstreamES，downstreamEE。如下图所示

![skipped_exon](https://i.loli.net/2021/01/18/YFJ2a7zG4mf3td6.png)

* 3. 13-16列  展示两组样品在IJC(inclusion junction) 和SJC(skipping junction counts)下的counts数，重复样本的结果以逗号分隔：列名分别为IJC_SAMPLE_1，SJC_SAMPLE_1，IJC_SAMPLE_2，SJC_SAMPLE_2。如下图

![exon](https://i.loli.net/2021/01/18/gkyW4SUwREcnB5F.png)

* 4. lncFormLen和SkipFormLen分别是inclusion form和skipping form的有效长度值,虽然有计算公式，还是要根据reads跨越时的具体情况来定。
* 5. IncLevel 可被认作为exon inclusion level(φ),是exon inclusion isoform在总(Exon inclusion isoform +exon skipping isoform 所占比例)
* 6. IncLevelDifference则是指两组样本IncLevel的差异，如果一组内多个样本，那么则是各自组的均值之间差值




references to [rMATS差异可变剪切分析](https://cloud.tencent.com/developer/article/1366294)
