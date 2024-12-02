# Faceoffs
NHL Player Faceoff Rankings 23/24
### Percentage of faceoffs won is a widely valued stat in ice hockey. Winning a faceoff allows a team to establish possession after a stoppage of play, and in the NHL the leaders in faceoff percentage often correlate closely with the ranks of the most reputable defensive forwards in the league - whether the faceoff itself is a crucial defensive act or not. While the percentage won is a useful summary, it fails to take into account the quality of opposition - if two players had the same skill but one faced much harder competition, that player would lose a higher percentage of their faceoffs than the other - and that players lose far more of their faceoffs on their so-called "off side" - with the forehand side of their stick facing toward the boards - such that many coaches will even have one left-handed forward take nearly all faceoffs on the left side of the ice, and another right-handed forward take nearly all faceoffs on the right side.

### We can take these variables into account and create an interpretable measure of a player's skill at faceoffs using a Bradley-Terry model: pairwise logistic regression where the probability of success is player 1's coefficient minus player 2's coefficient. We fit this model in Stan, a flexible probabilistic programming language for fitting Bayesian models. 
