# !!! IMPORTANT !!! Disable Antialias and Rounded line ends in plot window

if(ARG1 eq '')\
    print "usage: gnuplot [-p] [-e \"out='png'\"] -c ".ARG0." datafile [lowerBound] [upperBound]";\
    print "use -e \"out='png'\" for huge data sets. png can be replaced with whatever terminal your gnuplot supports";\
    exit

datafile=ARG1

# bounds for partial plotting for numbers >2^32 floating point notation is needed (e.g. 4.2)
lowerBound=ARG2
upperBound=ARG3

# use png for huge data sets
if(exists('out'))\
    set term out;\
    set output datafile[:strstrt(datafile,'.')].out

unset key

set xlabel 'CPU cycles'
set ylabel 'Thread id'

#set xtics 200000000

# set ranges
# required to fix offsets -> wtf
set auto fix
set offsets 0, 0, 0.5, 0.5
if(lowerBound ne '')\
    if(upperBound ne '')\
        set xrange [lowerBound:upperBound];\
    else\
        set xrange [lowerBound:]

set palette model RGB
set palette defined ( 0.0 '#000F89',\
                      0.5 '#000F89',\
                      0.5 '#F6B200',\
                      1.5 '#F6B200',\
                      1.5 '#84F600',\
                      2.5 '#84F600',\
                      2.5 '#FF6700',\
                      3.5 '#FF6700',\
                      3.5 '#00DCB7',\
                      4.5 '#00DCB7',\
                      4.5 '#FFF700',\
                      5.5 '#FFF700',\
                      5.5 '#870034',\
                      6.5 '#870034',\
                      6.5 '#7F01E0',\
                      7.0 '#7F01E0')
set cbrange [0:7]
set cbtics (\
            "P2P\nPreprocess\nSingleCell"  0,\
            "P2P\nc08Step\nBlock"          2,\
            "P2P\nPostprocess\nSingleCell" 1,\
            "M2L\nInitialize\nCell"        3,\
            "M2L\nInitialize\nSource"      4,\
            "M2L\nTranslation"             6,\
            "M2L\nFinalize\nCell"          5,\
           )

plot datafile \
     using \
        (lowerBound eq '' ?\
            column('Start') :\
            (column('Start') >= lowerBound ?\
                (upperBound eq '' ?\
                    column('Start') :\
                    (column('Stop') <= upperBound ?\
                        column('Start') :\
                        1/0\
                    )\
                ) :\
                1/0\
            )\
        ) :\
        'Thread':\
        'Start' :\
        'Stop'  :\
        (column('Thread')-.2):\
        (column('Thread')+.2):\
        'Type'\
     with boxxy \
     fillstyle solid .75 \
     palette
