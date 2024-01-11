Tube test analysis tutorial
================
Justin Varholick, PhD
January 11, 2024

## Background

The tube test is a commonly used behavioral assay to help determine the
dominance ranks of mice or rats kept in the laboratory. Briefly, a
narrow tube is placed horizontally on a table and cage-mates are placed
at opposite ends of the tube. The cage-mates will meet in the middle and
the animal that pushes the other out is considered the “winner” of that
trial. All cage-mates compete with one another multiple times in a
round-robin tournament to determine who is alpha, beta, gamma, epsilon,
etc. More details on the protocol can be found from the [Nature Methods
paper](https://doi.org/10.1038/s41596-018-0116-4), and I published an
[overview and
critique](https://www.researchgate.net/publication/352464789_Competitive_Exclusion)
in 2019.

<figure>
<img src="other%20resources/tube_test.gif"
alt="A gif of the tube test in action between two mice. A door is removed from the middle of the tube to start the interaction. The test is performed under red-light during the dark-cycle since mice are nocturnal." />
<figcaption aria-hidden="true">A gif of the tube test in action between
two mice. A door is removed from the middle of the tube to start the
interaction. The test is performed under red-light during the dark-cycle
since mice are nocturnal.</figcaption>
</figure>

Below is the general pipeline/method I used to determine the dominance
ranks of animals using the tube-test. I used this method in the
following publications: <https://doi.org/10.1038/s41598-018-24624-4> and
<https://doi.org/10.1038/s41598-019-49612-0>

You can download this R Markdown file and run the same code with your
own data.  
Then check out my other tutorials on graphing dominance data (forthcoming)

## Summary of data analytics used

This work uses the R packages, dplyr and tidyr to organize the data into
a usable form. The data is then run with the ‘compete’ package to
determine the dominance ranks. Custom graphs were made with ggplot to
visualize the changes in dominance rank across three weeks. The data is
then filtered and analyzed per behavioral test (exploration, zero-maze,
stress hormones, and body mass) with complementary plots.

## Determining the dominance ranks with the ‘compete’ package

### Packages

First, we will load the packages into R.

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
library(compete) #devtools::install_github('jalapic/compete') 
```

### Import data and clean

Let’s first import the data and look at the structure

``` r
#import data
tube_test_outcomes <- read.csv("data/tube_test_outcomes.csv", nrows=546)
head(tube_test_outcomes)
```

    ##   batch test_day experimenter cage_id opponent_1 opponent_2 t1 t2 t3 t4
    ## 1     1        1         Dean      41          1          5  2  2  2  1
    ## 2     1        1         Dean      41          2          3  1  1  1  2
    ## 3     1        1         Dean      41          4          5  1  1  1  1
    ## 4     1        1         Dean      41          1          3  2  1  1  2
    ## 5     1        1         Dean      41          2          4  2  2  2  2
    ## 6     1        1         Dean      41          3          5  1  1  2  2
    ##   avg_outcome winner start_side_opponent_1
    ## 1        1.75    0.0                     R
    ## 2        1.25    1.0                     R
    ## 3        1.00    1.0                     R
    ## 4        1.50    0.5                     L
    ## 5        2.00    0.0                     L
    ## 6        1.50    0.5                     R

The most important columns are the cage_id, opponent_1, opponent_2 and
the trails; t1, t2, t3, t4. The opponent columns refer to the individual
animal_id, who is competing within whom (e.g., animal 1 vs. animal 5).
The trial columns refer to which animal won the interaction, 1 refers to
opponent_1 and 2 refers to opponent_2. The other columns (e.g., batch,
test_day, experimenter) are less important, but can be used to test for
differences in ranks between experimenters, batches, days, etc.

We now need to restructure the data for analysis.

1.  First, I would like to concatenate test_day and cage_id so I can
    calculate a different dominance rank per test_day.

2.  Second, I would like the opponents to have clear animal_ids, so I
    need to concatenate cage_id and opponent_x.

3.  Third, I would like a long data set where the outcomes of all trials
    are in the same column.

``` r
#1.concatenate cage_id and test_day
tube_test_outcomes <- tube_test_outcomes %>%
  select(test_day,cage_id, opponent_1, opponent_2, t1, t2, t3, t4) %>%
  within(cage_day <- paste(cage_id, test_day, sep = "-d")) %>%
  select(opponent_1, opponent_2, t1, t2, t3, t4, cage_day, cage_id)
#2.rename opponent data to animal_id
tube_test_outcomes$opponent_1 <- paste(tube_test_outcomes$cage_id, tube_test_outcomes$opponent_1, sep="-")
tube_test_outcomes$opponent_2 <- paste(tube_test_outcomes$cage_id, tube_test_outcomes$opponent_2, sep="-")
#3.make a long dataset
long_data <- data.frame(tube_test_outcomes[1:2], tube_test_outcomes[7:8], stack(tube_test_outcomes[3:6]))
long_data <- long_data[,c("opponent_1", "opponent_2", "values", "ind", "cage_day", "cage_id")]
head(long_data)
```

    ##   opponent_1 opponent_2 values ind cage_day cage_id
    ## 1       41-1       41-5      2  t1    41-d1      41
    ## 2       41-2       41-3      1  t1    41-d1      41
    ## 3       41-4       41-5      1  t1    41-d1      41
    ## 4       41-1       41-3      2  t1    41-d1      41
    ## 5       41-2       41-4      2  t1    41-d1      41
    ## 6       41-3       41-5      1  t1    41-d1      41

### Convert data tables into lists

Now the data is in a long format where values refer to the outcome of
the tube test trial and ind refers to the corresponding trial.

To determine the dominance ranks for each cage_id or cage_day, we need
to divide the data into a list. Then we can apply a function to the list
using lapply and append the results to a single data frame using rbind.

``` r
#split up long_data by day_id
day_list <- split(long_data, as.factor(long_data$cage_day))

#splitup long_dat by Cage (compound all 3 trials)
cage_list <- split(long_data, as.factor(long_data$cage_id))
```

### David’s Score

Now, we can calculate the David’s score using the lists and lapply with
rbind.

``` r
#make a custom function that formats the table into a win-loss matrix and then runs the ds function.
get_ds <- function(x){
  matrix <- get_wl_matrix(x, ties="keep")
  david_scores <- as.data.frame(ds(matrix, norm=FALSE, type = "D"))
}

#get david scores for each item in the list
ds_per_week <- do.call(rbind, lapply(day_list, get_ds))

#tidy data
ds_per_week <- ds_per_week %>% 
  rename(ds = `ds(matrix, norm = FALSE, type = "D")`) %>% 
  mutate(cage_day = rownames(.)) %>% 
  separate_wider_delim(cage_day, delim = ".", names = c("cage_day", "animal_id")) %>% 
  separate_wider_delim(cage_day, delim = "-", names = c("cage_id", "day")) %>% 
  select(cage_id, animal_id, day, ds)

write.csv(ds_per_week, "data/ds_per_week.csv")

head(ds_per_week)
```

    ## # A tibble: 6 × 4
    ##   cage_id animal_id day       ds
    ##   <chr>   <chr>     <chr>  <dbl>
    ## 1 110     110-1     d1     7    
    ## 2 110     110-2     d1     3    
    ## 3 110     110-3     d1    -4    
    ## 4 110     110-4     d1    -5.5  
    ## 5 110     110-5     d1    -0.500
    ## 6 110     110-1     d2     6

We can use this same method for the cage_list

``` r
#get david scores for each item in the list
ds_per_cage <- do.call(rbind, lapply(cage_list, get_ds))

#tidy data
ds_per_cage <- ds_per_cage %>% 
  rename(ds = `ds(matrix, norm = FALSE, type = "D")`) %>% 
  mutate(cage_animal_id = rownames(.)) %>% 
  separate_wider_delim(cage_animal_id, delim = ".", names = c("cage_id", "animal_id")) %>% 
  select(cage_id, animal_id, ds)

write.csv(ds_per_week, "data/ds_per_cage.csv")

head(ds_per_cage)
```

    ## # A tibble: 6 × 3
    ##   cage_id animal_id    ds
    ##   <chr>   <chr>     <dbl>
    ## 1 38      38-2      -2.31
    ## 2 38      38-3       1.23
    ## 3 38      38-4      -1.85
    ## 4 38      38-5       2.92
    ## 5 41      41-1       2.31
    ## 6 41      41-2      -3.85

### Dominance Rank

The same method can be used to determine the dominance rank, and any
other metric

``` r
get_MinRank <- function(x){
  matrix <- get_wl_matrix(x, ties="keep")
  wins <- rowSums(matrix)
  total_wins <- as.data.frame(wins)
  total_wins$ID <- rownames(total_wins)
  total_wins <- mutate(total_wins, rank = rank(-wins, ties.method="min"))
  total_wins <- select(total_wins, ID, wins, rank)
  return(total_wins)
}

#ranks per week
ranks_per_week <- do.call(rbind, lapply(day_list, get_MinRank))
#tidy
ranks_per_week <- ranks_per_week %>% 
  mutate(cage_day = rownames(.)) %>% 
  separate_wider_delim(cage_day, delim = ".", names = c("cage_day", "animal_id")) %>% 
  separate_wider_delim(cage_day, delim = "-", names = c("cage_id", "day")) %>% 
  select(cage_id, animal_id, day, wins, rank)
write.csv(ranks_per_week, "data/ranks_per_week.csv")

#ranks per cage
ranks_per_cage <- do.call(rbind, lapply(cage_list, get_MinRank))
#tidy
ranks_per_cage <- ranks_per_cage %>% 
  mutate(tempid = rownames(.)) %>% 
  separate_wider_delim(tempid, delim = ".", names = c("cage_id", "animal_id")) %>% 
  select(cage_id, animal_id, wins, rank)
write.csv(ranks_per_cage, "data/ranks_per_cage.csv")
```

### Landau’s h

Then use the same general method to calculate the Landau’s score (the
linearity of the dominance hierarchy)

``` r
Landau <- function(x){
  matrix <- get_wl_matrix(x, ties="remove")
  devries(matrix)
}

#each day
landau_per_week <- do.call(rbind, lapply(day_list, Landau))
landau_per_week <- as.data.frame(landau_per_week) %>% 
  mutate(tempid = rownames(.)) %>% 
  separate_wider_delim(tempid, delim = "-", names = c("cage_id", "day"))
landau_per_week$`h-modified` <- as.numeric(landau_per_week$`h-modified`)
landau_per_week$`p-value` <- as.numeric(landau_per_week$`p-value`)
write.csv(as.data.frame(landau_per_week), "data/landau_per_week.csv")
#each cage
landau_per_cage <- do.call(rbind, lapply(cage_list, Landau))
landau_per_cage <- as.data.frame(landau_per_cage) %>% 
  mutate(cage_id = rownames(.))
landau_per_cage$`h-modified` <- as.numeric(landau_per_cage$`h-modified`)
landau_per_cage$`p-value` <- as.numeric(landau_per_cage$`p-value`)
write.csv(landau_per_cage, "data/landau_per_cage.csv")
```

### Stability

We can also restructure the data to determine whether individual animals
were stable or unstable across test_days

``` r
#add animal_id to ranks_per_week
animal_inventory <- read.csv("data/animal_inventory.csv")
sex_id <- select(animal_inventory, animal_id, sex)
ranks_per_week <- left_join(ranks_per_week, sex_id, by = "animal_id")
#add the cage_rank to ranks_per_week
cage_rank <- ranks_per_cage %>% 
  select(animal_id, rank) %>% 
  rename(cage_rank = rank)
ranks_per_week <- left_join(ranks_per_week, cage_rank, by = "animal_id")

stability <- ranks_per_week %>% 
  select(cage_id, animal_id, sex, cage_rank, day, rank) %>% 
  pivot_wider(names_from = day, values_from = rank) %>% 
  mutate(stability = case_when(
    d1 == d2 & d2 == d3 ~ "stable",
    d1 != d2 ~ "unstable",
    d1 != d3 ~ "unstable",
    d2 != d3 ~ "unstable"
  )) %>% 
  filter(stability == "stable")

write.csv(stability, "data/stability.csv")

head(stability)
```

    ## # A tibble: 6 × 8
    ##   cage_id animal_id sex    cage_rank    d1    d2    d3 stability
    ##   <chr>   <chr>     <chr>      <int> <int> <int> <int> <chr>    
    ## 1 110     110-1     Male           1     1     1     1 stable   
    ## 2 113     113-1     Male           5     5     5     5 stable   
    ## 3 113     113-2     Male           2     3     3     3 stable   
    ## 4 113     113-3     Male           1     1     1     1 stable   
    ## 5 122     122-3     Female         4     4     4     4 stable   
    ## 6 122     122-5     Female         3     3     3     3 stable
