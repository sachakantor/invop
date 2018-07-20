import math

deliveryPoints = [11, 47, 44, 25, 10, 26, 26, 54, 18, 51, 20, 105, 7, 16, 43, 100, 50, 21, 11, 19, 14, 10, 11]
spirit = [34, 411, 82, 157, 5, 183, 14, 215, 102, 21, 54, 0, 6, 96, 118, 112, 535, 8, 53, 28, 69, 65, 27]
oilFirstRegion = [9, 13, 14, 17, 18, 19, 23, 21]
oilSecondRegion = [9, 11, 17, 18, 18, 17, 22, 24, 36, 43]
oilThirdRegion = [6, 15, 15, 25, 39]
growthCategoryA = [1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0]
growthCategoryB = [0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1]

def printDeliveryPointsConstraints():
    printSummationConstraints(deliveryPoints)

def printSpiritMarketControlConstraints():
    printSummationConstraints(spirit)

def printOilMarketRegionControlConstraints():
    printSummationConstraints(oilFirstRegion)
    printSummationConstraints(oilSecondRegion, 8)
    printSummationConstraints(oilThirdRegion, 18)

def printGrowthMarketConstraints():
    printSummationConstraints(growthCategoryA)
    printSummationConstraints(growthCategoryB)

def printDeliveryPointsDeviationConstraints():
    printStandardDeviationConstraints(deliveryPoints)

def printSpiritMarketControlDeviationConstraints():
    printStandardDeviationConstraints(spirit)

def printOilMarketRegionControlDeviationConstraints():
    printStandardDeviationConstraints(oilFirstRegion)
    printStandardDeviationConstraints(oilSecondRegion, 8)
    printStandardDeviationConstraints(oilThirdRegion, 18)

def printGrowthMarketDeviationConstraints():
    printStandardDeviationConstraints(growthCategoryA)
    printStandardDeviationConstraints(growthCategoryA)

def printStandardDeviationConstraints(coefficients, retailerNumberOffset = 0):
    summation = ""
    for i in range(0, len(coefficients)):
        summation += str(coefficients[i]) + " M_" + str(i + 1 + retailerNumberOffset)
        if i != len(coefficients) - 1:
            summation += " + "
    coefficientsSum = sum(coefficients)
    print(summation + " + " + str(coefficientsSum) + " U => " + str(coefficientsSum * 0.4))
    print(summation + " - " + str(coefficientsSum) + " U <= " + str(coefficientsSum * 0.4))

def printSummationConstraints(coefficients, retailerNumberOffset = 0):
    summation = ""
    for i in range(0, len(coefficients)):
        summation += str(coefficients[i]) + " M_" + str(i + 1 + retailerNumberOffset)
        if i != len(coefficients) - 1:
            summation += " + "
        coefficientsSum = sum(coefficients)
    print(summation + " >= " + str(coefficientsSum * 0.35))
    print(summation + " <= " + str(coefficientsSum * 0.45))

def printBinaryVariables():
    for retailer in range(1, 24):
        print("M_" + str(retailer))

print("MINIMIZE")
print("U")
print("ST.")
printDeliveryPointsConstraints()
printSpiritMarketControlConstraints()
printOilMarketRegionControlConstraints()
printGrowthMarketConstraints()
printDeliveryPointsDeviationConstraints()
printSpiritMarketControlDeviationConstraints()
printOilMarketRegionControlDeviationConstraints()
printGrowthMarketDeviationConstraints()
print("BOUNDS")
print("BINARY")
printBinaryVariables()
print("END")