#!/bin/gnuplot -p
# !!! IMPORTANT !!! Disable Antialias and Rounded line ends in plot window

if(ARG1 eq '')\
    print "usage: gnuplot [-p] [-e \"out='png'\"] -c ".ARG0." datafile [lowerBound] [upperBound]";\
    print "use -e \"out='png'\" for huge data sets. png can be replaced with whatever terminal your gnuplot supports";\
    print "Large bounds need to be entered as floats (42.0)";\
    exit

datafile=ARG1

# bounds for partial plotting
lowerBound=ARG2
upperBound=ARG3

# use png for huge data sets
if(exists('out'))\
    set term out;\
    set output datafile[0:strstrt(datafile,'.')].out

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
set palette defined ( 0.0  '#011BE1',\
                      0.5  '#011BE1',\
                      0.5  '#001089',\
                      1.5  '#001089',\
                      1.5  '#000949',\
                      2.5  '#000949',\
                      2.5  '#6F0083',\
                      3.5  '#6F0083',\
                      3.5  '#C72900',\
                      4.5  '#C72900',\
                      4.5  '#FFB800',\
                      5.5  '#FFB800',\
                      5.5  '#F6B200',\
                      6.5  '#F6B200',\
                      6.5  '#9C7100',\
                      7.5  '#9C7100',\
                      7.5  '#C79000',\
                      8.5  '#C79000',\
                      8.5  '#C79000',\
                      9.5  '#C79000',\
                      9.5  '#ABC200',\
                      10.5 '#ABC200',\
                      10.5 '#008E3E',\
                      11.0 '#008E3E')
set cbrange [0:11]
set cbtics (\
            "P2P\nPreprocess"               0,\
            "P2P\nc08Step"                  2,\
            "P2P\nPostprocess"              1,\
            "P2M"                           3,\
            "M2M"                           4,\
            "M2L\nInit Cell"                5,\
            "M2L\nInit Source"              6,\
            "M2L\nCompleteCell"             8,\
            "M2L\nPair2Way"                 9,\
            "M2L\nFinalize"                 7,\
            "L2L"                          10,\
            "L2P"                          11,\
           )
# set cbtics (\
            # "P2P\nPreprocess\nSingleCell"   0,\
            # "P2P\nc08Step\nBlock"           2,\
            # "P2P\nPostprocess\nSingleCell"  1,\
            # "P2M\nCompleteCell"             3,\
            # "M2M\nCompleteCell"             4,\
            # "M2L\nInitialize\nCell"         5,\
            # "M2L\nInitialize\nSource"       6,\
            # "M2L\nTranslation"              8,\
            # "M2L\nPair2Way"                 9,\
            # "M2L\nFinalizeCell"             7,\
            # "L2L\nCompleteCell"            10,\
            # "L2P\nCompleteCell"            11,\
           # )

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

