def likelyhood_ratio(baseline, observed):
    if observed < 0:
        raise ValueError("The number of observed events is negative")
    if baseline < 0:
        raise ValueError("The number of simulated events is negative")
    if baseline == 0 and observed > 0:
        return 1.0
    if observed > baseline:
        return np.power((observed / baseline), observed) * np.exp(baseline - observed)
    return 1.0
