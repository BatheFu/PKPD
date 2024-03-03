# A two Roll Dice with R

roll <- function(){
    die <- 1:6;
    dice <- sample(die,size = 2, replace = T);
    sum(dice)
}

roll()

roll2 <- function(bones){
    dice <- sample(bones, size = 2, replace = T);
    sum(dice)
}

roll2(bones = 1:6)

rolls <- replicate(10000, roll()) %>% 
qplot(binwidth = 1)

roll3 <- function() {
  die <- 1:6;
  dice <- sample(die, size = 2, replace = T,
                 prob = c(1/8,1/8,1/8,1/8,1/8,3/8));
  sum(dice);
}

replicate(10000, roll3()) %>% 
    qplot(binwidth = 1)

