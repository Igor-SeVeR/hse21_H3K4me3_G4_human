
# hse21_H3K4me3_G4_human
## Финальный проект по майнору «Биоинформатика»
### Егоров Игорь Сергеевич, группа МБ-2

### Часть 1

#### Выбранные структуры

| Организм | Структура ДНК | Гистоновая метка | Тип клеток | Метка 1 | Метка 2 |
| -------- | ------------- | ---------------- | ---------- | ------- | ------- |
| Human (hg19) | G4_seq_Li_K | H3K4me3 | H1 | [ENCFF668YOE](https://www.encodeproject.org/files/ENCFF668YOE/) | [ENCFF254ACI](https://www.encodeproject.org/files/ENCFF254ACI/) |

#### Подготовка репозитория

Для начала создадим свой [github](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human) репозиторий для данного проекта. Далее склонируем данный репозиторий на сервер и создадим 3 основных папки:
[images](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/tree/main/images) - здесь расположены все картинки (графики, скриншоты и т. д.) полученные в результате работы;
[src](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/tree/main/src) - здесь расположены все скрипты, написанные во время выполнения проекта;
[data](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/tree/main/data) - здесь расположены все файлы, с которыми велась работа на протяжении всего проекта.
Запушим все изменения на github, таким образом получим красиво оформленный, структурированный репозиторий.

#### Анализ пиков гистоновой метки

##### Подготовка файлов

Переходим в папку data и скачиваем туда 2 .bed файла ChIP-seq экспериментов из ENCODE:

```bash
wget https://www.encodeproject.org/files/ENCFF254ACI/@@download/ENCFF254ACI.bed.gz
wget https://www.encodeproject.org/files/ENCFF668YOE/@@download/ENCFF668YOE.bed.gz
```

Далее оставим а данных файлах только первые 5 столбцов:

```bash
zcat ENCFF254ACI.bed.gz  |  cut -f1-5 > H3K4me3_H1.ENCFF254ACI.hg38.bed
zcat ENCFF668YOE.bed.gz  |  cut -f1-5 > H3K4me3_H1.ENCFF668YOE.hg38.bed
```

Мы скачали данные файлы в формате hg38, нам необхоимо перегнать данные файлы в формат hg18. Для этого мы скачиваем файл с сопоставлением координат одного формата и другого:

```bash
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz
```

После этого перегоняем наши файлы из одного формата в другой:

```bash
liftOver   H3K4me3_H1.ENCFF254ACI.hg38.bed   hg38ToHg19.over.chain.gz   H3K4me3_H1.ENCFF254ACI.hg19.bed   H3K4me3_H1.ENCFF254ACI.unmapped.bed
liftOver   H3K4me3_H1.ENCFF668YOE.hg38.bed   hg38ToHg19.over.chain.gz   H3K4me3_H1.ENCFF668YOE.hg19.bed   H3K4me3_H1.ENCFF668YOE.unmapped.bed
```

Проверим корректность преобразования файлов. Для этого посмотрим на количество строк в исходных файлах, преобразованных файлах, а также в unmapped файлах:

```bash
isegorov@laboratory01:~/ngs/final_project$ wc -l H3K4me3_H1.ENCFF254ACI.hg38.bed
26245 H3K4me3_H1.ENCFF254ACI.hg38.bed
isegorov@laboratory01:~/ngs/final_project$ wc -l H3K4me3_H1.ENCFF668YOE.hg38.bed
34171 H3K4me3_H1.ENCFF668YOE.hg38.bed
isegorov@laboratory01:~/ngs/final_project$ wc -l H3K4me3_H1.ENCFF254ACI.hg19.bed
26165 H3K4me3_H1.ENCFF254ACI.hg19.bed
isegorov@laboratory01:~/ngs/final_project$ wc -l H3K4me3_H1.ENCFF668YOE.hg19.bed
34105 H3K4me3_H1.ENCFF668YOE.hg19.bed
isegorov@laboratory01:~/ngs/final_project$ wc -l H3K4me3_H1.ENCFF254ACI.unmapped.bed
160 H3K4me3_H1.ENCFF254ACI.unmapped.bed
isegorov@laboratory01:~/ngs/final_project$ wc -l H3K4me3_H1.ENCFF668YOE.unmapped.bed
132 H3K4me3_H1.ENCFF668YOE.unmapped.bed
```

Как видим, числа не сходятся: для
ENCFF668YOE: 34171 - 34105 != 132;
ENCFF254ACI: 26245 - 26165 != 160.
Это достаточно легко объяснить - для этого изучим строение файлов umapped:


```bash
isegorov@laboratory01:~/ngs/hse21_H3K4me3_G4_human/data$ head H3K4me3_H1.ENCFF254ACI.unmapped.bed
#Partially deleted in new
chr1    120197261       120197624       Peak_13752      71
#Partially deleted in new
chr1    121052087       121052476       Peak_13754      71
#Split in new
chr1    145424453       145425166       Peak_3896       270
#Deleted in new
chr1    145426765       145427089       Peak_13760      71
#Partially deleted in new
chr1    145427250       145430609       Peak_2612       386
```

Как видим каждая запись занимает две строчки - причина, по которой запись не была сматчена и сама запись.
Имея это знание, получаем иные числа:
ENCFF668YOE: 34171 - 34105 = 66;
ENCFF254ACI: 26245 - 26165 = 80.
Отлично, всё сошлось.
Остаётся только запушить все файлы в github репозиторий.

##### Построение гистограмм длин участков

Клонируем репозитрой на локальный компьютер, где и будет производиться работа по созданию данного проекта.
Напишем [скрипт](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/src/len_hist.R) для построения гистограмм длин участков для каждого эксперимента. Далее с помощью него построим графики для уже полученных четырёх файлов:

![len_hist.H3K4me3_H1.ENCFF254ACI.hg19](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/len_hist.H3K4me3_H1.ENCFF254ACI.hg19-1.png)

![len_hist.H3K4me3_H1.ENCFF254ACI.hg38](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/len_hist.H3K4me3_H1.ENCFF254ACI.hg38-1.png)

![len_hist.H3K4me3_H1.ENCFF668YOE.hg19](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/len_hist.H3K4me3_H1.ENCFF668YOE.hg19-1.png)

![len_hist.H3K4me3_H1.ENCFF668YOE.hg38](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/len_hist.H3K4me3_H1.ENCFF668YOE.hg38-1.png)

Количество пиков можно видеть на графиках.

##### Фильтрация ChIP-seq

Теперь выкинем из ChIP-seq пиков слишком длинные пики (outliers). Для этого отсортируем [H3K4me3_H1.ENCFF254ACI.hg19](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.ENCFF254ACI.hg19.bed), [H3K4me3_H1.ENCFF668YOE.hg19](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.ENCFF668YOE.hg19.bed) и посмотрим на самые длинные пики:

```bash
H3K4me3_H1.ENCFF668YOE.hg19
   chrom     start       end       name score  len
1  chr12 113904197 113913711 Peak_12479   172 9514
2  chr20  47893547  47901069   Peak_184   647 7522
3   chr9  23819631  23826983  Peak_7732   264 7352
4  chr17  47072021  47079367   Peak_618   555 7346
5  chr20  17945679  17952993   Peak_339   603 7314
6  chr11  65265233  65272419   Peak_584   560 7186
7  chr12  95938999  95946147  Peak_4227   361 7148
8  chr14  69256662  69263711   Peak_965   520 7049
9   chr3  32858059  32864977  Peak_1214   501 6918
10  chr2  64066889  64073672   Peak_877   530 6783
```

Видим, что в [H3K4me3_H1.ENCFF668YOE.hg19](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.ENCFF668YOE.hg19.bed) есть один аутлейнер - рид, длина которого 9514, его нужно удалить.


```bash
Longest reads
H3K4me3_H1.ENCFF254ACI.hg19
   chrom     start       end     name score    len
1   chr7  27154493  27257439   Peak_1  1000 102946
2   chr2 176972442 177047064  Peak_38  1000  74622
3  chr17  46655470  46727684  Peak_27  1000  72214
4  chr12  54378950  54449016  Peak_68  1000  70066
5  chr12  54336250  54378676  Peak_67  1000  42426
6  chr17  46792856  46833921  Peak_13  1000  41065
7   chr1 119520088 119551591 Peak_111  1000  31503
8   chr2 172942634 172974000  Peak_95  1000  31366
9  chr13  28528189  28558578  Peak_11  1000  30389
10  chr8  11550415  11580674 Peak_195  1000  30259
```

А вот в [H3K4me3_H1.ENCFF254ACI.hg19](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.ENCFF254ACI.hg19.bed) аутлейнеров немного больше - 4, их длины 50000+.

Для удаления аутлейнеров воспользуемся [скриптом](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/src/len_filter.R). В итоге получим два файла: [H3K4me3_H1.ENCFF254ACI.hg19.filtered.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.ENCFF254ACI.hg19.filtered.bed), [H3K4me3_H1.ENCFF668YOE.hg19.filtered.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.ENCFF668YOE.hg19.filtered.bed).

Теперь построим гистограммы длин участков для этих двух файлов:

![len_hist.H3K4me3_H1.ENCFF254ACI.hg19.filtered](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/len_hist.H3K4me3_H1.ENCFF254ACI.hg19.filtered-1.png)

![len_hist.H3K4me3_H1.ENCFF668YOE.hg19.filtered](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/len_hist.H3K4me3_H1.ENCFF668YOE.hg19.filtered-1.png)

Количество пиков можно видеть на графиках.

#####  Расположение пиков гистоновой метки относительно аннотированных генов

Для наглядности построим график типа пай-чарт с помощью R-библиотеки ChIPseeker и библиотеки с аннотацией генов (разметкой).
Напишем [скрипт](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/src/chip_seeker.R) позволяющий нам по bed файлу получить такую диаграмму.

Далее прогоним на нём оба наших отфильтрованных файла. В итоге получаем:
для [H3K4me3_H1.ENCFF254ACI.hg19.filtered.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.ENCFF254ACI.hg19.filtered.bed):
![chip_seeker.H3K4me3_H1.ENCFF254ACI.hg19.filtered.plotAnnoPie](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/chip_seeker.H3K4me3_H1.ENCFF254ACI.hg19.filtered.plotAnnoPie-1.png)
для [H3K4me3_H1.ENCFF668YOE.hg19.filtered.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.ENCFF668YOE.hg19.filtered.bed):
![chip_seeker.H3K4me3_H1.ENCFF668YOE.hg19.filtered.plotAnnoPie](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/chip_seeker.H3K4me3_H1.ENCFF668YOE.hg19.filtered.plotAnnoPie-1.png)

Также для каждого отфильтрованного файла построим дополнительно график, который будет показывать, где пики находятся на хромосомах. Это очень информативный и полезный график, однако, из-за огромного количества пиков он крайне плохо интерпретируется.
Для [H3K4me3_H1.ENCFF254ACI.hg19.filtered.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.ENCFF254ACI.hg19.filtered.bed):
![chip_seeker.H3K4me3_H1.ENCFF254ACI.hg19.filtered.covplot](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/chip_seeker.H3K4me3_H1.ENCFF254ACI.hg19.filtered.covplot-1.png)
Для [H3K4me3_H1.ENCFF668YOE.hg19.filtered.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.ENCFF668YOE.hg19.filtered.bed):
![chip_seeker.H3K4me3_H1.ENCFF668YOE.hg19.filtered.covplot](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/chip_seeker.H3K4me3_H1.ENCFF668YOE.hg19.filtered.covplot-1.png)

##### Объединение двух наборов отфильтрованных ChIP-seq пиков

Для выполнения данного шага, нам понадобится перенести данные на сервер. Для это пушим все изменения с локальной машины на гитхаб, а на сервере просто подтягиваем изменения репозитория.
Теперь при помощи утилиты bedtools merge объединим два наших отфильтрованных файла в один:

```bash
cat  *.filtered.bed  |   sort -k1,1 -k2,2n   |   bedtools merge   >  H3K4me3_H1.merged.hg19.bed
```

Отлично, мы получили файл [H3K4me3_H1.merged.hg19.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.merged.hg19.bed). Давайте посмотрим на количество строк в нём:

```bash
isegorov@laboratory01:~/ngs/hse21_H3K4me3_G4_human/data$ wc -l H3K4me3_H1.merged.hg19.bed
48006 H3K4me3_H1.merged.hg19.bed
```

Как видим, в данном файле строк меньше, чем суммарно в двух файлах, из которых он образован. Это связано с наличием пересечений, которые bedtools merge удаляет.

Построим для [H3K4me3_H1.merged.hg19.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.merged.hg19.bed) две аналогичные предыдущему пункту диаграммы:

![chip_seeker.H3K4me3_H1.merged.hg19.plotAnnoPie](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/chip_seeker.H3K4me3_H1.merged.hg19.plotAnnoPie-1.png)

![chip_seeker.H3K4me3_H1.merged.hg19.covplot](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/chip_seeker.H3K4me3_H1.merged.hg19.covplot-1.png)

По графикам можно заметить, что данные на них ожидаемы (для пайчарта - нечто среднее между графиками двух файлов, из которых образован [H3K4me3_H1.merged.hg19.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.merged.hg19.bed), для второго же графика - их объединение).

##### Визуализация исходных наборов ChIP-seq пиков

Кульминацией первой части проекта является загрузка всех трёх файлов: [H3K4me3_H1.merged.hg19.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.merged.hg19.bed), [H3K4me3_H1.ENCFF668YOE.hg19.filtered.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.ENCFF668YOE.hg19.filtered.bed), [H3K4me3_H1.ENCFF254ACI.hg19.filtered.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.ENCFF254ACI.hg19.filtered.bed)
на [genome.ucsc.edu](genome.ucsc.edu), где мы визуализируем наши наборы ChIP-seq пиков, а также проверим корректность merge.
Загрузка данных на [genome.ucsc.edu](genome.ucsc.edu) осуществляется при помощи следующих команд:

```bash
track visibility=dense name="ENCFF254ACI"  description="H3K4me3_H1.ENCFF254ACI.hg19.filtered.bed"
https://raw.githubusercontent.com/Igor-SeVeR/hse21_H3K4me3_G4_human/main/data/H3K4me3_H1.ENCFF254ACI.hg19.filtered.bed
track visibility=dense name="ENCFF668YOE"  description="H3K4me3_H1.ENCFF668YOE.hg19.filtered.bed"
https://raw.githubusercontent.com/Igor-SeVeR/hse21_H3K4me3_G4_human/main/data/H3K4me3_H1.ENCFF668YOE.hg19.filtered.bed
track visibility=dense name="ChIP_merge"  color=50,50,200   description="H3K4me3_H1.merged.hg19.bed"
https://raw.githubusercontent.com/Igor-SeVeR/hse21_H3K4me3_G4_human/main/data/H3K4me3_H1.merged.hg19.bed
```

Все три файла загружаем по очереди во вкладку [My Data-> Custom Tracks](http://genome.ucsc.edu/cgi-bin/hgCustom). Впрочем, это было более, чем детально, рассмотрено на семинаре, так что расписывать данный пункт не вижу особого смысла.
Далее создаём [сессию](http://genome.ucsc.edu/s/isegorov/H3K4me3_H1), в которой и сохраняем нашу визуализацию трёх файлов. При желании, уважаемый читатель может перейти в неё и поиграться.

Для доказательства корректности merge приведу скриншот:

![Checking_merging_correctness_genome_ucsc_edu.png](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/Checking_merging_correctness_genome_ucsc_edu.png)

Как видим из скриншота: объединение полностью покрывает оба набора ChIP-seq пиков. При большом желании можно походить по хромосоме в сессии и убедиться в этом окончательно.

### Часть 2

#### Анализ участков вторичной структуры ДНК

##### Получение вторичной структцры ДНК

Необходимо скачать соответствующую нашему эксперименту вторичную структуру ДНК. В моём слуае это структура, соответствующая G4_seq_Li_K. Нахоим её по [ссылке из таблицы c распределением](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3003539).
Там нас ждут два файла [Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed) и [Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed).
Нам нужно их скачать.

```bash
isegorov@laboratory01:~/ngs/hse21_H3K4me3_G4_human/data$ wc -l Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed
216321 Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed
isegorov@laboratory01:~/ngs/hse21_H3K4me3_G4_human/data$ wc -l Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed
217951 Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed
```

Однако, на данный момент в них слишком много ненужной нам информации. Обрежем ненужные столбцы:

```bash
zcat GSM3003539_Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed.gz | cut -f1-5 > GSM3003539_Homo_minus.bed
zcat GSM3003539_Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed.gz | cut -f1-5 > GSM3003539_Homo_plus.bed
```

Отлично! Всё что нам осталось - это объединить два данных файла:

```bash
cat GSM3003539_Homo_minus.bed GSM3003539_Homo_plus.bed | sort -k1,1 -k2,2n | bedtools merge > GSM3003539_Homo.bed
```

Получаем файл [GSM3003539_Homo.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/GSM3003539_Homo.bed), который и представляет из себя вторичную структуру ДНК.
Теперь поанализируем данные файлы. Посмотрим на длину файлов перед объединением и после:

```bash
isegorov@laboratory01:~/ngs/hse21_H3K4me3_G4_human/data$ wc -l Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed
216321 Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed
isegorov@laboratory01:~/ngs/hse21_H3K4me3_G4_human/data$ wc -l Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed
217951 Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed
isegorov@laboratory01:~/ngs/hse21_H3K4me3_G4_human/data$ wc -l GSM3003539_Homo.bed
428624 GSM3003539_Homo.bed
```

Как видим, длина объединенённого файла вновь получилось меньше, чем сумма длин тех файлов, из которых он был образован. Это снова связано с удалением дубликатов.

(Здесь и далее я не буду расписывать простейшие действия, связанные с поддержанием актуального github репозитория. Будем считать, что все полученные новые файлы сразу отправляются в него)

##### Распределение длин участков вторичной структуры ДНК

С помощью уже знакомого нам [скрипта](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/src/len_hist.R) строим гистограмму:

![len_hist.GSM3003539_Homo](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/len_hist.GSM3003539_Homo-1.png)

Количество пиков указано на графике.

##### Расположение участков структуры ДНК относительно аннотированных генов

С данной операцией мы тоже уже хорошо знакомы, так что запускаем [скрипт](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/src/chip_seeker.R) на [GSM3003539_Homo.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/GSM3003539_Homo.bed) и получаем пайчарт:

![chip_seeker.GSM3003539_Homo.plotAnnoPie](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/chip_seeker.GSM3003539_Homo.plotAnnoPie-1.png)

##### Визуализация вторичной структуры ДНК

Как и в первой части, кульминацией второй части проекта является подгрузка файла со вторичной структурой ДНК на [genome.ucsc.edu](genome.ucsc.edu), где она визуализируется.

```bash
track visibility=dense name="GSM3003539_Homo"  color=100,100,100   description="GSM3003539_Homo.bed"
https://raw.githubusercontent.com/Igor-SeVeR/hse21_H3K4me3_G4_human/main/data/GSM3003539_Homo.bed
```

В итоге получаем новую [сессию](http://genome.ucsc.edu/s/isegorov/H3K4me3_H1_GSM3003539_Homo), где уже можно поизучать и местонахождение вторичной структуры ДНК.

![Added_secondary_dna_genome_ucsc_edu](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/Added_secondary_dna_genome_ucsc_edu.png)

Фрагмент полученной визуализации.

### Часть 3

#### Анализ пересечений гистоновой метки и структуры ДНК

##### Построение пересечения гистоновых меток с структурой ДНК

Найти пересечение гистновых меток и вторичной структуры ДНК нам позволит bedtools intersect:

```bash
bedtools intersect -a GSM3003539_Homo.bed -b H3K4me3_H1.merged.hg19.bed > H3K4me3_H1.intersect_with_GSM3003539_Homo.bed
```

Получаем файл [H3K4me3_H1.intersect_with_GSM3003539_Homo.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.intersect_with_GSM3003539_Homo.bed). Сразу изучим его длину.

```bash
39812 H3K4me3_H1.intersect_with_GSM3003539_Homo.bed
48006 H3K4me3_H1.merged.hg19.bed
```

Данное число пересечений является однозначно великолепным. Оно составляет свыше 0,8 всех пиков эксперимента, что достаточно много.

##### Распределение длин пересечений

Вновь возвращаемся к уже хорошо знакомому нам [скрипту](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/src/len_hist.R), с помощью которого строим гистограмму длин пересечений:

![len_hist.H3K4me3_H1.intersect_with_GSM3003539_Homo](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/len_hist.H3K4me3_H1.intersect_with_GSM3003539_Homo-1.png)

Количество пиков указано на графике.

Также дополнительно построим пайчарт, показывающий нам расположение пиков пересечения относительно аннотированных генов:

![chip_seeker.H3K4me3_H1.intersect_with_GSM3003539_Homo.plotAnnoPie](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/chip_seeker.H3K4me3_H1.intersect_with_GSM3003539_Homo.plotAnnoPie-1.png)

##### Визуализация в геномном браузере

Подгружаем наш [H3K4me3_H1.intersect_with_GSM3003539_Homo.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.intersect_with_GSM3003539_Homo.bed) в [genome.ucsc.edu](genome.ucsc.edu):

```bash
track visibility=dense name="intersect_with_GSM3003539_Homo"  color=255,0,0  description="H3K4me3_H1.intersect_with_GSM3003539_Homo.bed"
https://raw.githubusercontent.com/Igor-SeVeR/hse21_H3K4me3_G4_human/main/data/H3K4me3_H1.intersect_with_GSM3003539_Homo.bed
```

Снова формируем новую [сессию](http://genome.ucsc.edu/s/isegorov/H3K4me3_GSM3003539_Intersection), которая является для проекта финальной.
Поищем в ней одно-два места, где есть пересечение между гистоновой меткой и стр-рой ДНК. Их несложно найти:
первое место:
![Intersection_secondary_struct_with_histone_marks_genome_ucsc_edu_first](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/Intersection_secondary_struct_with_histone_marks_genome_ucsc_edu_first.png)
координаты пересечения = chr1:78,957,334-78,957,438;
второе место:
![Intersection_secondary_struct_with_histone_marks_genome_ucsc_edu_second](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/Intersection_secondary_struct_with_histone_marks_genome_ucsc_edu_second.png)
координаты пересечения = chr1:75,598,259-75,598,363.

##### Ассоциирование полученных пересечений с ближайшими генами

С помощью R-библиотеки ChIPpeakAnno напишем [скрипт](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/src/ChIPpeakAnno.R), который позволит нам ассоциировать полученные пересечения с ближайшими к ним генами.
После выполнения скрипта получим два файла: [H3K4me3_H1.intersect_with_GSM3003539_Homo.genes.txt](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.intersect_with_GSM3003539_Homo.genes.txt), [H3K4me3_H1.intersect_with_GSM3003539_Homo.genes_uniq.txt](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.intersect_with_GSM3003539_Homo.genes_uniq.txt).
Изучим их строение.

```bash
isegorov@laboratory01:~/ngs/hse21_H3K4me3_G4_human/data$ head H3K4me3_H1.intersect_with_GSM3003539_Homo.genes_uniq.txt
LINC01128
SAMD11
KLHL17
NOC2L
PLEKHN1
PERM1
HES4
ISG15
AGRN
C1orf159
```

Делаем закономерный вывод - файл [H3K4me3_H1.intersect_with_GSM3003539_Homo.genes_uniq.txt](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.intersect_with_GSM3003539_Homo.genes_uniq.txt) содержит перечисление всех уникальных генов, с которыми были ассоциированы пересечения. Давайте выясним, сколько их всего:

```bash
isegorov@laboratory01:~/ngs/hse21_H3K4me3_G4_human/data$ wc -l  H3K4me3_H1.intersect_with_GSM3003539_Homo.genes_uniq.txt
9171 H3K4me3_H1.intersect_with_GSM3003539_Homo.genes_uniq.txt
```

Получаем 9171 уникальный ген.

Теперь изучим второй файл:

```bash
isegorov@laboratory01:~/ngs/hse21_H3K4me3_G4_human/data$ head H3K4me3_H1.intersect_with_GSM3003539_Homo.genes.txt
seqnames        start   end     width   strand  peak    feature feature.ranges.start    feature.ranges.end      feature.ranges.width    feature.strand  distance        insideFeature       distanceToSite  symbol
chr1    762044  762071  28      *       X00001  643837  762971  794826  31856   +       899     upstream        899     LINC01128
chr1    860270  860419  150     *       X00011  148398  860530  879961  19432   +       110     upstream        110     SAMD11
chr1    894601  894710  110     *       X00033  339451  895967  901099  5133    +       1256    upstream        1256    KLHL17
chr1    894601  894710  110     *       X00033  26155   879583  894679  15097   -       0       overlapStart    0       NOC2L
chr1    895441  896055  615     *       X00034  339451  895967  901099  5133    +       0       overlapStart    0       KLHL17
chr1    895441  896055  615     *       X00034  26155   879583  894679  15097   -       761     upstream        761     NOC2L
chr1    900929  901048  120     *       X00035  84069   901877  910484  8608    +       828     upstream        828     PLEKHN1
chr1    901189  901408  220     *       X00036  84069   901877  910484  8608    +       468     upstream        468     PLEKHN1
chr1    919360  919448  89      *       X00042  84808   910579  917473  6895    -       1886    upstream        1886    PERM1
```

Снова делаем вывод: файл [H3K4me3_H1.intersect_with_GSM3003539_Homo.genes.txt](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.intersect_with_GSM3003539_Homo.genes.txt) содежит в себе первой строкой - описание формата, а далее список всех пересечений, каждому из которых ассоциирован определённый ген.
Посмотрим на его длину:

```bash
isegorov@laboratory01:~/ngs/hse21_H3K4me3_G4_human/data$ wc -l  H3K4me3_H1.intersect_with_GSM3003539_Homo.genes.txt
15984 H3K4me3_H1.intersect_with_GSM3003539_Homo.genes.txt
```

То есть количество пиков с ассоциированными им генами равно 15984 - 1(строка с описанием формата) = 15983

##### GO-анализ для полученных уникальных генов

GO - генная онтология. Это анализ генов, который позвляет предположить, к какой "системе" относится данный набор генов (понятие "система" станет понятно после вывода результатов GO-анализа). Это вероятностный тест, использущий в своей основе критерий Фишера.
Для проведения данного анализа воспользумся [ресурсом](http://pantherdb.org/). Как проводить сам эксперимент на данном ресурсе было детально описано в [инструкции выполнения проекта](https://docs.google.com/document/d/1wtbNQ0ZN3ruHHIg5-PgxuGyX2TiIoJSD2GkzR-AyeCI/edit#), так что данные шаги я описывать не буду.
Загружаем наши уникальные гены, ставим эксперимент. В итоге получаем следующий результат:

![GO_H3K4me3_H1.intersect_with_GSM3003539_Homo.genes_uniq_result](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/GO_H3K4me3_H1.intersect_with_GSM3003539_Homo.genes_uniq_result.png)

У нас немного unmapped или multiple mapped генов, что не может не радовать. Также получаем [таблицу](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/GO_H3K4me3_H1.intersect_with_GSM3003539_Homo.genes_uniq_analysis.txt) с предсказанными категориями. Расмотрим самые значимые из них:

|GO biological process complete | # | # | expected | Fold Enrichment | +/- | raw P value | FDR |
| -------- | ------------- | ---------------- | ---------- | ------- | ------- | ------- | ------- |
detection of chemical stimulus involved in sensory perception|486|0|200.68|< 0.01|-|2.54E-74|4.02E-70
detection of chemical stimulus involved in sensory perception of smell|441|0|182.10|< 0.01|-|3.43E-67|2.71E-63
detection of chemical stimulus|522|13|215.54|.06|-|7.40E-61|3.90E-57
sensory perception of chemical stimulus|543|26|224.21|.12|-|9.61E-52|3.80E-48
sensory perception of smell|468|17|193.24|.09|-|2.22E-49|7.01E-46
detection of stimulus involved in sensory perception|557|40|229.99|.17|-|1.45E-43|3.82E-40
regulation of cellular metabolic process|6138|3188|2534.48|1.26|+|7.82E-37|1.77E-33
regulation of primary metabolic process|5925|3081|2446.53|1.26|+|2.01E-35|3.97E-32
positive regulation of cellular process|5741|2996|2370.55|1.26|+|5.22E-35|9.16E-32
regulation of nitrogen compound metabolic process|5740|2991|2370.14|1.26|+|1.70E-34|2.68E-31
regulation of metabolic process|6874|3480|2838.38|1.23|+|5.49E-34|7.89E-31
detection of stimulus|721|101|297.71|.34|-|1.38E-31|1.81E-28

Данная таблица была получена путём сортировки значения FDR в результах эксперимента на том же [ресурсе](http://pantherdb.org/)
