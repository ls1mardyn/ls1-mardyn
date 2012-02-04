import csv
import os

with open('domainInfos.csv', 'rb') as domainInfoFile:
    with open('domainHistogram.csv', 'wb') as domainHistogramFile:
        domainInfoReader = csv.reader(domainInfoFile)
        domainHistogramWriter = csv.writer(domainHistogramFile)

        maxSize = 0

        domainInfos = []
        for row in domainInfoReader:
            config = os.path.splitext(os.path.basename(row.pop(0)))[0];
            numCells = len(row)
            histogram = {}
            for cellCount in row:
                cellCount = int(cellCount)
                if cellCount in histogram:
                    histogram[cellCount] = histogram[cellCount] + 1
                else:
                    histogram[cellCount] = 1

            sortedCounts = histogram.keys()
            sortedCounts.sort()
            occurences = [histogram[count] for count in sortedCounts]
            maxSize = max([maxSize, len(sortedCounts)])

            plotY = sum([[1,y,y,1] for y in sortedCounts],[0])
            """plotX = [0]
            for i in range(len(sortedCounts)):
                x1 = plotX[-1]
                x2 = x1 + sortedCounts[i]*occurences[i]
                plotX = plotX + [x1,x1,x2,x2]"""
                
            plotX = reduce( lambda l, d: l + [l[-1], l[-1], l[-1]+d, l[-1]+d], occurences, [0] )
            maxWidth = float(plotX[-1])
            plotX = [ x/maxWidth for x in plotX ]

            domainInfos.append( {
                'config':config,
                'cellMoleculeCount':sortedCounts,
                'frequency':occurences,
                'size':len(sortedCounts),
                'plotX':plotX,
                'plotY':plotY
                } )

        """domainHistogramWriter.writerow( ['plottable data:'] )
        for domainInfo in domainInfos:
            filenameColumn = [domainInfo['config']]
            domainHistogramWriter.writerow( filenameColumn + domainInfo['plotX'] )
            domainHistogramWriter.writerow( filenameColumn + domainInfo['plotY'] )"""

        for domainInfo in domainInfos:
            domainHistogramWriter.writerow( [domainInfo['config'] + '.plotX'] + domainInfo['plotX'] )
            domainHistogramWriter.writerow( [domainInfo['config'] + '.plotY'] + domainInfo['plotY'] )
            domainHistogramWriter.writerow( [domainInfo['config'] + '.cellMoleculeCount'] + domainInfo['cellMoleculeCount'] )
            domainHistogramWriter.writerow( [domainInfo['config'] + '.frequency'] + domainInfo['frequency'] )

        """row = []
        for domainInfo in domainInfos:
            row.extend( [ domainInfo['config'] + '.cellMoleculeCount', domainInfo['config'] + '.frequency' ] )
        domainHistogramWriter.writerow( row )
        for i in range(maxSize):
            row = []
            for domainInfo in domainInfos:
                if i < domainInfo['size']:
                    row.extend( [domainInfo['cellMoleculeCount'][i], domainInfo['frequency'][i]] )
                else:
                    row.extend( [None, None] )

            domainHistogramWriter.writerow( row )"""
                    

            
