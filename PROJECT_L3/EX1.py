
import math 

def basic_tm(S):
    Tm = 4 * (S.count('G') + S.count('C')) + 2 * (S.count('A') + S.count('T'))
    return Tm

def advanced_tm(S):
    na = 0.001
    length = len(S)
    gc_count = S.count('G') + S.count('C')
    gc_perc = (gc_count / length) * 100
    log_na = math.log10(na)
    tm = 81.5 + 16.6 * log_na + 0.41 * gc_perc - (600 / length)
    return round(tm, 2)

S = 'TCCAGACGACTA'

print('The basic formula for Tm [Tm = 4 * (G + C) + 2 * (A + T)] has the following result:', basic_tm(S), '°C')

print('\nThe advanced formula for Tm [Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) – 600/length] has the following result:', advanced_tm(S), '°C')
