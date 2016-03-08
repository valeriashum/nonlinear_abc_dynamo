# Start of statistical function definitions
def mean(list):
    return sum(list) / len(list)

def variance(list):
    squaredList = [x*x for x in list]
    return sum(squaredList)/len(list) - (sum(list)/len(list))**2

def covariance(list1, list2):
    if len(list1) != len(list2):
        return 0
    else:
        return sum([x*y for x,y in zip(list1, list2)]) / len(list1) -  sum(list1)*sum(list2)/(len(list1)**2)
# End of statistical function definitions
