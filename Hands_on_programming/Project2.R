# Playing Cards

# use R to assemble a deck of 52 playing cards.
rm(list = ls())

dice <- 1:6;
is.vector(dice) # TRUE

# An attribute is a piece of information that you can attach to an atomic vector (or any R object).
# The attribute won’t affect any of the values in the object, 
# and it will not appear when you display your object.
# You can think of an attribute as “metadata”;

# The most common attributes to give an atomic vector
# are names, dimensions (dim), and classes.

names(dice)  <- c("one", "two", "three", "four", "five", "six")

attributes(dice)
#the names won’t affect the actual values of the vector

hand1 <- c("ace", "king", "queen", "jack", "ten", "spades", "spades", 
           "spades", "spades", "spades")

# Three ways for matrix
matrix(hand1, nrow = 5)
matrix(hand1, ncol = 2)
dim(hand1) <- c(5, 2)

hand2 <- c("ace", "spades", "king", "spades", "queen", "spades", "jack", 
           "spades", "ten", "spades")
# Fill the matrix by row, instead of default column 
matrix(hand2, nrow = 5, byrow = TRUE)
matrix(hand2, ncol = 2, byrow = TRUE)

Sys.time()

deck <- data.frame(
    face = c("king", "queen", "jack", "ten", "nine", "eight", "seven", "six",
             "five", "four", "three", "two", "ace", "king", "queen", "jack", "ten", 
             "nine", "eight", "seven", "six", "five", "four", "three", "two", "ace", 
             "king", "queen", "jack", "ten", "nine", "eight", "seven", "six", "five", 
             "four", "three", "two", "ace", "king", "queen", "jack", "ten", "nine", 
             "eight", "seven", "six", "five", "four", "three", "two", "ace"),  
    suit = c("spades", "spades", "spades", "spades", "spades", "spades", 
             "spades", "spades", "spades", "spades", "spades", "spades", "spades", 
             "clubs", "clubs", "clubs", "clubs", "clubs", "clubs", "clubs", "clubs", 
             "clubs", "clubs", "clubs", "clubs", "clubs", "diamonds", "diamonds", 
             "diamonds", "diamonds", "diamonds", "diamonds", "diamonds", "diamonds", 
             "diamonds", "diamonds", "diamonds", "diamonds", "diamonds", "hearts", 
             "hearts", "hearts", "hearts", "hearts", "hearts", "hearts", "hearts", 
             "hearts", "hearts", "hearts", "hearts", "hearts"), 
    value = c(13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 13, 12, 11, 10, 9, 8, 
              7, 6, 5, 4, 3, 2, 1, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 13, 12, 11, 
              10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
)


deck[1,c(1,2,3)]

