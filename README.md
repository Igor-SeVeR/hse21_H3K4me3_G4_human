
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

####  Расположение пиков гистоновой метки относительно аннотированных генов

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

#### Объединение двух наборов отфильтрованных ChIP-seq пиков

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

По графикам можно заметить, что данные на них ожидаемы (нечто среднее между графиками двух файлов, из которых образован [H3K4me3_H1.merged.hg19.bed](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.merged.hg19.bed))

#### Визуализация исходных наборов ChIP-seq пиков

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

Все 3 файла загружаем по очереди во вкладку [My Data-> Custom Tracks](http://genome.ucsc.edu/cgi-bin/hgCustom). Впрочем, это было более, чем детально рассмотрено на семинаре, так что расписывать данный пункт не вижу особого смысла.
Далее создаём [сессию](http://genome.ucsc.edu/s/isegorov/H3K4me3_H1), в которой и сохраняем нашу визуализацию трёх файлов. При желании, уважаемый читатель может перейти в неё и поиграться.

Для доказательства корректности merge приведу скриншот:

![Checking_merging_correctness_genome_ucsc_edu.png](https://github.com/Igor-SeVeR/hse21_H3K4me3_G4_human/blob/main/images/Checking_merging_correctness_genome_ucsc_edu.png)

Как видим из скриншота: объединение полностью покрывает оба набора ChIP-seq пиков. При большом желании можно походить по хромосоме в сессии и убедиться в этом окончательно.
