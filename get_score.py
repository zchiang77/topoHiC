import time
import logging
from collections import Counter
import pandas as pd
import numpy as np


def get_score(exp_name, directory, score_threshold, bin_range, res=None, loop_lim=None):
    logging.basicConfig(filename=exp_name + '.log', level=logging.INFO,
                        format='%(asctime)s - %(name)s - %(levelname)s:%(message)s')
    logger = logging.getLogger(__name__)

    logger.info("Step4: Computing loop enrichment score.")
    logger.info("Score > {} will be recognized as TADs.".format(score_threshold))

    start_time = time.time()

    listLoops = []
    lengthLoops = []
    i = 0
    with open(directory + 'scratch.txt', 'r') as file:
        for line in file:
            i += 1
            line = line.strip('\n')
            line = line.split(',')
            line = line[:-1]
            line = [int(x) for x in line]
            if loop_lim is None:
                listLoops.append(line)
                lengthLoops.append(max(line) - min(line))
            else:
                if (max(line) - min(line)) < loop_lim / res:
                    listLoops.append(line)
                    lengthLoops.append(max(line) - min(line))

    infoLoopsCounter = Counter(lengthLoops)
    infoLoops = pd.DataFrame(infoLoopsCounter.items(), columns=['length', 'number'])
    start = bin_range[0]
    end = bin_range[1]

    score = pd.DataFrame(columns=['bin_id', 'loop_enrichment', 'weight'])
    for key in range(start, end + 1):
        enrichment = 0
        size = 0
        for loop in listLoops:
            # minAnchor = min(loop)
            # maxAnchor = max(loop)
            # if (key == minAnchor) or (key == maxAnchor):
            #     enrichment += 1
            #     size += maxAnchor - minAnchor
            if key in loop:
                enrichment += 1
                size += max(loop) - min(loop)
        if size == 0:
            weight = 0
        else:
            weight = np.mean(size)
        score.loc[len(score)] = [key, enrichment, weight]

    logger.info("Loop enrichment score computation completed.")

    score['bin_id'] = score['bin_id'].astype(int)
    score.to_csv(directory + 'score.csv', header=True, index=False)
    infoLoops.to_csv(directory + 'loopInfo.csv', header=True, index=False)

    end_time = time.time()
    logger.info("Step4 completed! Time taken(s): {:.2}.".format(end_time - start_time))

    return score

