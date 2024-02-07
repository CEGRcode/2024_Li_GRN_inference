import numpy as np
import pandas as pd

def Combine_redundant_attractors(input_path, output_path):
    RNAseq_data = open('input_path, 'r')
    StableStatesCollector = {}
    TranscriptionPofiles = []
    for line in RNAseq_data:
        line_split = line.split()
        if line_split[0] == 'ABF1':
            continue
        else:
            for i in range(0, len(line_split)):
                line_split[i] = float(line_split[i])
        TranscriptionPofiles.append(line_split[1:])
        StableStatesCollector[str(line_split[1:])] = 1

    TranscriptionPofiles = np.array(TranscriptionPofiles)
    TranscriptionPofileMax = np.max(TranscriptionPofiles, axis=0)

    TransMaxOver1 = []
    for each in TranscriptionPofileMax:
        TransMaxOver1.append(1/each)

    VectorCollector = []
    FreqCollector = []
    for keys in StableStatesCollector:
        if StableStatesCollector[keys] >= 0: # Threshold for the Stable States frequency
            VectorCollector.append(keys)
            FreqCollector.append(StableStatesCollector[keys])
        else:
            pass
    for i in range(0, len(VectorCollector)):
        VectorCollector[i] = VectorCollector[i][:-1]
        VectorCollector[i] = VectorCollector[i][1:]
    MatrixCollector = np.random.randint(2,size=(1, len(VectorCollector), len(VectorCollector[0].split())))[0]
    for i in range(0, len(VectorCollector)):
        TempVec = VectorCollector[i].split()
        for j in range(0, len(TempVec)):
            MatrixCollector[i][j] = float(TempVec[j][0:-1])

    MyCounter = 0
    StableStateThreshold = 0.3
    #print(len(MatrixCollector))
    while MyCounter < MatrixCollector.shape[0] - 1:
        j = MyCounter + 1
        if MyCounter == MatrixCollector.shape[0] - 1 - 1:
            #print(np.dot(abs(MatrixCollector[MyCounter] - MatrixCollector[j]), TransMaxOver1))
            if np.dot(abs(MatrixCollector[MyCounter] - MatrixCollector[j]), TransMaxOver1) < StableStateThreshold:
                print(MatrixCollector[MyCounter], '\n', MatrixCollector[j], '\n', np.dot(abs(MatrixCollector[MyCounter] - MatrixCollector[j]), TransMaxOver1), '\n\n')
                TempCombineVector = (MatrixCollector[MyCounter] + MatrixCollector[j]) * 0.5
                MatrixCollector[MyCounter] = TempCombineVector
                MatrixCollector = np.delete(MatrixCollector, j, axis = 0)
                TempCombineFreq = FreqCollector[MyCounter] + FreqCollector[j]
                FreqCollector[MyCounter] = TempCombineFreq
                FreqCollector.pop(j)
            else:
                pass
            break
        else:
            while True:
                if np.dot(abs(MatrixCollector[MyCounter] - MatrixCollector[j]), TransMaxOver1) < StableStateThreshold:
                    print(MatrixCollector[MyCounter], '\n', MatrixCollector[j], '\n', np.dot(abs(MatrixCollector[MyCounter] - MatrixCollector[j]), TransMaxOver1), '\n\n')
                    TempCombineVector = (MatrixCollector[MyCounter] + MatrixCollector[j]) * 0.5
                    MatrixCollector[MyCounter] = TempCombineVector
                    MatrixCollector = np.delete(MatrixCollector, j, axis = 0)
                    TempCombineFreq = FreqCollector[MyCounter] + FreqCollector[j]
                    FreqCollector[MyCounter] = TempCombineFreq
                    FreqCollector.pop(j)
                    break
                else:
                    j = j + 1
                if j == MatrixCollector.shape[0]:
                    MyCounter = MyCounter + 1
                    break
                else:
                    pass

    #print(StableStatesCollector)

    for i in range(0, len(MatrixCollector)):
        print(str(MatrixCollector[i]))
    print(len(MatrixCollector))

    outfile = open(output_path, 'a')
    for each in MatrixCollector:
        if len(np.argwhere(each == 0)) == 0:
               outfile.write('-1\t')
        else:
               outfile.write(str(np.argwhere(each == 0)[0][0])+'\t')
        for values in each:
               outfile.write(str(values)+'\t')
        outfile.write('\n')
    outfile.close()
    return
