import time
import logging
import pandas as pd
from scipy.stats import shapiro


def my_concat(set1, set2):
    if len(set1) == 0 and len(set2) == 0:
        return set1
    elif len(set1) == 0 or len(set2) == 0:
        if len(set1) == 0:
            return set2
        else:
            return set1
    else:
        return pd.concat([set1, set2], axis=0)

def my_shapiro(all_start, all_end, domain_file, score_file):
    p_list = []
    if len(domain_file) > 1:
        for i in range(0, len(domain_file)):
            if i == 0:
                sw_start = all_start
                sw_end = domain_file.loc[1]['start'] - 1
            elif i == len(domain_file) - 1:
                sw_start = domain_file.loc[i - 1]['end'] + 1
                sw_end = all_end
            else:
                sw_start = domain_file.loc[i - 1]['end'] + 1
                sw_end = domain_file.loc[i + 1]['start'] - 1
            sw_start = min([sw_start, domain_file.loc[i]['start'] - 4])
            sw_end = max([sw_end, domain_file.loc[i]['end'] + 4])
            test_seq = score_file[((score_file['bin_id'] >= sw_start) & (score_file['bin_id'] <= sw_end))][
                'loop_enrichment']
            stat, p = shapiro(test_seq)
            p_list.append(p)
    else:
        sw_start = all_start
        sw_end = all_end
        test_seq = score_file[((score_file['bin_id'] >= sw_start) & (score_file['bin_id'] <= sw_end))][
            'loop_enrichment']
        stat, p = shapiro(test_seq)
        p_list.append(p)
    domain_file['pvalue'] = p_list
    return domain_file

def iter_sep_domain(sub_score, cand_start, cand_end, thresh, lim):
    sep_domain = pd.DataFrame(columns=['start', 'end', 'length', 'threshold'])
    cand_start = int(cand_start)
    cand_end = int(cand_end)
    search_sub_bin = cand_start
    while search_sub_bin <= cand_end:
        sub_enr = sub_score[sub_score['bin_id'] == search_sub_bin]['loop_enrichment'].tolist()[0]
        if sub_enr > thresh:
            subdomain_start = search_sub_bin
            for search_sub_end in range(subdomain_start, cand_end + 1):
                search_sub_enr = sub_score[sub_score['bin_id'] == search_sub_end]['loop_enrichment'].tolist()[0]
                if search_sub_enr > thresh:
                    continue
                else:
                    break
            subdomain_end = search_sub_end - 1
            if subdomain_end - subdomain_start >= lim:
                sep_domain.loc[len(sep_domain)] = [subdomain_start, subdomain_end,
                                                   subdomain_end - subdomain_start, thresh]
            search_sub_bin = search_sub_end + 1
        else:
            search_sub_bin += 1
            continue
    return sep_domain

def get_domain(directory, score_threshold, bin_range, test_sig, max_itr, res, lim, meg):

    start = bin_range[0]
    end = bin_range[1]
    limit = int(lim / res)

    domain = pd.DataFrame(columns=['start', 'end', 'length', 'pvalue', 'threshold'])
    if meg:
        score = pd.read_csv(directory + '_score.csv')
    else:
        score = pd.read_csv(directory + '/score.csv')

    tmp_threshold = score_threshold

    tmp_domain = iter_sep_domain(sub_score=score, cand_start=start, cand_end=end,
                                 thresh=tmp_threshold, lim=limit)
    tmp_domain = my_shapiro(all_start=start, all_end=end, domain_file=tmp_domain, score_file=score)
    domain = my_concat(domain, tmp_domain[tmp_domain['pvalue'] >= test_sig])
    domain.reset_index(drop=True, inplace=True)

    itr = 0
    while len(tmp_domain) != 0:

        itr += 1
        if (itr > max_itr) and (max_itr >= 0):
            break

        if len(tmp_domain[tmp_domain['pvalue'] < test_sig]) == 0:
            break

        next_domain = tmp_domain[tmp_domain['pvalue'] < test_sig]
        next_domain.reset_index(drop=True, inplace=True)

        tmp_domain = pd.DataFrame()
        tmp_threshold += 1
        for j in range(0, len(next_domain)):
            next_start = next_domain.loc[j]['start']
            next_end = next_domain.loc[j]['end']
            next_subdomain = iter_sep_domain(sub_score=score, cand_start=next_start,
                                             cand_end=next_end, thresh=tmp_threshold,
                                             lim=limit)
            if len(next_subdomain) != 0:
                next_subdomain = my_shapiro(next_start, next_end, next_subdomain, score)
                if len(next_subdomain[next_subdomain['pvalue'] >= test_sig]) > 0:
                    domain = my_concat(domain, next_subdomain[next_subdomain['pvalue'] >= test_sig])
                    domain.sort_values(by='start', inplace=True)
                    domain.reset_index(drop=True, inplace=True)
                tmp_domain = my_concat(tmp_domain, next_subdomain[next_subdomain['pvalue'] < test_sig])
        tmp_domain.reset_index(drop=True, inplace=True)

    if max_itr == 0:
        domain = tmp_domain
    else:
        domain = my_concat(domain, tmp_domain)

    domain.sort_values(by='start', inplace=True)
    domain.reset_index(drop=True, inplace=True)
    domain[['start', 'end']] = domain[['start', 'end']].astype(int)
    domain['pvalue'] = domain['pvalue'].round(3)

    if meg:
        domain.to_csv(directory + '_domain.csv', header=True, index=False)
    else:
        domain.to_csv(directory + '/domain.csv', header=True, index=False)

    return domain




