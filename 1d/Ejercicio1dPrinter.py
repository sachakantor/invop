import math

typeToAmount = [("1", 12), ("2", 10), ("3", 5)]
demandRanges = [(0, 6), (6, 9), (9, 15), (15, 18), (18, 24)]
demands = [15000, 30000, 25000, 40000, 27000]
typeToLevels = {"1": ("850", "2000"), "2": ("1250", "1750"), "3": ("1500", "4000")}
typeToMinimumCost = {"1": 1000, "2": 2600, "3": 3000}
typeToCostPerHour = {"1": 2, "2": 1.3, "3": 3}
typeToStartCost = {"1": 2000, "2": 1000, "3": 500}

def printDemandConstraints():
    for index, (start, end) in enumerate(demandRanges):
        constraint = []
        hourlyRange = range(start, end)
        enumerateDemandVariables(
            hourlyRange, lambda demandVar, t, i, h: constraint.append(demandVar + " + ")
        )
        for h in hourlyRange:
            constraint.append("900 A_" + str(h) + " + ")
        for h in hourlyRange:
            constraint.append("1400 B_" + str(h))
            if h != end - 1:
                constraint.append(" + ")
        constraint.append(" >= " + str(demands[index]))
        print(''.join(constraint))

def printDemandGenerationLeftHandConstraints():
    enumerateDemandVariables(
        range(0, 24), lambda demandVar, t, i, h: print(typeToLevels[t][0] + " Y_" + t + i + h + " - " + demandVar + " <= 0")
    )

def printDemandGenerationRightHandConstraints():
    enumerateDemandVariables(
        range(0, 24), lambda demandVar, t, i, h: print(demandVar + " - " + typeToLevels[t][1] + " Y_" + t + i + h + " <= 0")
    )

def printHydroGenerationLeftHandConstraints():
    enumerateHydroVariables(
        range(0, 24), lambda hydroVar, t, i, h: print(typeToLevels[t][0] + " P_" + t + i + h + " - " + hydroVar + " <= 0")
    )

def printHydroGenerationRightHandConstraints():
    enumerateHydroVariables(
        range(0, 24), lambda hydroVar, t, i, h: print(hydroVar + " - " + typeToLevels[t][1] + " P_" + t + i + h + " <= 0")
    )

def enumerateDemandVariables(hourlyRange, receiver):
    for t, amount in typeToAmount:
        for i in range(1, amount + 1):
            for h in hourlyRange:
                receiver("W_" + t + str(i) + str(h), t, str(i), str(h))


def enumerateHydroVariables(hourlyRange, receiver):
    for t, amount in typeToAmount:
        for i in range(1, amount + 1):
            for h in hourlyRange:
                receiver("V_" + t + str(i) + str(h), t, str(i), str(h))

def printWaterDepthConstraints():
    constraint = ""
    for h in range(0, 24):
        constraint += "0.31 A_" + str(h) + " + 0.47 B_" + str(h)
        if h != 23:
            constraint += " + "
    constraint += " - "
    for h in range(0, 24):
        constraint += "P_" + str(h)
        if h != 23:
            constraint += " - "
    constraint += " = 0"
    print(constraint)

    printWaterDepthLeftHandConstraints()
    printWaterDepthRightHandConstraints()

def printWaterDepthLeftHandConstraints():
    for k in range(1, 23):
        constraint = ""
        for h in range(0, k + 1):
            constraint += "0.31 A_" + str(h) + " + 0.47 B_" + str(h)
            if h != k:
                constraint += " + "
        constraint += " - "
        for h in range(0, k + 1):
            constraint += "P_" + str(h)
            if h != k:
                constraint += " - "
        constraint += " <= 1"
        print(constraint)

def printWaterDepthRightHandConstraints():
    for k in range(1, 23):
        constraint = "- "
        for h in range(0, k + 1):
            constraint += "0.31 A_" + str(h) + " - 0.47 B_" + str(h)
            if h != k:
                constraint += " - "
        constraint += " + "
        for h in range(0, k + 1):
            constraint += "P_" + str(h)
            if h != k:
                constraint += " + "
        constraint += " <= 4"
        print(constraint)

def printFifteenPercentDemandIncreaseConstraints():
    for index, (start, end) in enumerate(demandRanges):
        hourlyRange = range(start, end)
        for h in hourlyRange:
            constraint = []
            constraint.append("- 900 A_" + str(h))
            constraint.append(" - 1400 B_" + str(h))
            enumerateDemandVariables(
                range(h, h + 1), 
                lambda demandVar, t, i, h: constraint.append(" + " + typeToLevels[t][1] + " Y_" + t + i + h + " - " + demandVar)
            )
            enumerateDemandVariables(
                range(h, h + 1),
                lambda demandVar, t, i, h: constraint.append(" + " + typeToLevels[t][1] + " P_" + t + i + h)
            )
            demandPercentage = math.ceil(0.15 * demands[index])
            constraint.append(" >= " + str(demandPercentage - 2300))
            print(''.join(constraint))

def printExclusiveGeneratorDedication():
    enumerateDemandVariables(
        range(0, 24), 
        lambda demandVar, t, i, h: print("P_" + t + i + h + " + Y_" + t + i + h + " < 2")
    )

def printHydroGenerationConstraints():
    for h in range(0, 24):
        constraint = []
        constraint.append("3000 P_" + str(h))
        enumerateHydroVariables(
            range(h, h + 1),
            lambda hydroVar, t, i, h: constraint.append(" - " + hydroVar)
        )
        constraint.append(" <= 0")
        print(''.join(constraint))

def printGeneratorsStartIfWereStopped():
    enumerateDemandVariables(
        range(1, 24),
        lambda demandVar, t, i, h: print("S_" + t + i + h + " + Y_" + t + i + str(int(h) - 1) + " < 2")
    )
    for h in range(1, 24):
        print("SA_" + str(h) + " + A_" + str(h - 1) + " < 2")
    for h in range(1, 24):
        print("SB_" + str(h) + " + B_" + str(h - 1) + " < 2")

def printGeneratorsUsageCoherency():
    enumerateDemandVariables(
        range(1, 24),
        lambda demandVar, t, i, h: print("Y_" + t + i + h + " - Y_" + t + i + str(int(h) - 1) + " - S_" + t + i + h + " <= 0")
    )
    enumerateDemandVariables(
        range(1, 24),
        lambda demandVar, t, i, h: print("P_" + t + i + h + " - P_" + t + i + str(int(h) - 1) + " - S_" + t + i + h + " <= 0")
    )

def printBinaryVariables():
    enumerateDemandVariables(
        range(0, 24),
        lambda var, t, i, h: print("Y_" + t + i + h)
    )
    enumerateDemandVariables(
        range(0, 24),
        lambda var, t, i, h: print("P_" + t + i + h)
    )
    enumerateDemandVariables(
        range(0, 24),
        lambda var, t, i, h: print("S_" + t + i + h)
    )
    dailyHourlyRange = range(0, 24)
    for h in dailyHourlyRange:
        print("A_" + str(h))
    for h in dailyHourlyRange:
        print("B_" + str(h))
    for h in dailyHourlyRange:
        print("SA_" + str(h))
    for h in dailyHourlyRange:
        print("SB_" + str(h))
    for h in dailyHourlyRange:
        print("P_" + str(h))

def printObjectiveFunction():
    def minCost(t): return typeToMinimumCost[t] - (2 * 850)
    dailyHourlyRange = range(0, 24)
    objectiveFunction = []
    for h in dailyHourlyRange:
        currentRange = range(h, h + 1)
        enumerateDemandVariables(
            currentRange,
            lambda demandVar, t, i, h: objectiveFunction.append(str(typeToCostPerHour[t]) + " " + demandVar + (" - " if minCost(t) < 0 else " + ") + str(abs(minCost(t))) + " Y_" + t + i + h + " + ")
        )
        enumerateHydroVariables(
            currentRange,
            lambda hydroVar, t, i, h: objectiveFunction.append(str(typeToCostPerHour[t]) + " " + hydroVar + (" - " if minCost(t) < 0 else " + ") + str(abs(minCost(t))) + " P_" + t + i + h + " + ")
        )
        enumerateDemandVariables(
            currentRange,
            lambda demandVar, t, i, h: objectiveFunction.append(str(typeToStartCost[t]) + " S_" + t + i + h + " + ")
        )
    for h in dailyHourlyRange:
        objectiveFunction.append("90 A_" + str(h) + " + ")
    for h in dailyHourlyRange:
        objectiveFunction.append("150 B_" + str(h) + " + ")
    for h in dailyHourlyRange:
        objectiveFunction.append("1500 SA_" + str(h) + " + ")
    for h in dailyHourlyRange:
        objectiveFunction.append("1200 SB_" + str(h))
        if (h != 23):
            objectiveFunction.append(" + ")
    print(''.join(objectiveFunction))

print("MINIMIZE")
printObjectiveFunction()
print("ST.")
printDemandConstraints()
printDemandGenerationLeftHandConstraints()
printDemandGenerationRightHandConstraints()
printHydroGenerationLeftHandConstraints()
printHydroGenerationRightHandConstraints()
printWaterDepthConstraints()
printFifteenPercentDemandIncreaseConstraints()
printExclusiveGeneratorDedication()
printGeneratorsStartIfWereStopped()
printHydroGenerationConstraints()
printGeneratorsUsageCoherency()
print("BOUNDS")
print("BINARY")
printBinaryVariables()
print("END")