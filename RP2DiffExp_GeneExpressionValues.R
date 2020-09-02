MaleExpression <- getEAWP(ymal)
FemaleExpression <- getEAWP(yfem)
PreExpression <- getEAWP(ypre)
PostExpression <- getEAWP(ypost)

write.csv(MaleExpression, 'MaleExpression.csv')
write.csv(FemaleExpression, 'FemaleExpression.csv')
write.csv(PreExpression, "PreExpression.csv")
write.csv(PostExpression, 'PostExpression.csv')
