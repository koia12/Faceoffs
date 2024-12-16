library(tidyverse)
library(rstan)
mc.cores <- parallel::detectCores()
#get rid of All-Star game (filtering by game id since it includes non-NHL players and would
#cause problems for inference to have players not connected to the network of player matchups)

faceoffs2 <- All_Events %>% filter(Event == "FAC" & as.numeric(str_remove(game_id, pattern = ("^[0-9]{4}"))) < 100000) %>% 
  select(c(player_1, player_2, homeTeamDefendingSide, homeTeam, awayTeam, xCoord,
           yCoord, zoneCode, away_number, home_number)) %>% mutate(player_1 = unlist(player_1),
                                                                   player_2 = unlist(player_2))


#we need to scrape whether player shoots left or right building off previous scraping
#get rid of All-Star game participants
boxscore_players_distinct <- boxscore_players %>% distinct(playerId, .keep_all = TRUE) %>% 
  filter(as.numeric(str_remove(game_id, pattern = ("^[0-9]{4}"))) < 100000)
boxscore_players_distinct$url <- str_c("https://api-web.nhle.com/v1/player/",
                                      boxscore_players_distinct$playerId,"/landing")


for (i in 1:nrow(boxscore_players_distinct)){
  boxscore_players_distinct$shootsCatches[i] <- jsonlite::fromJSON(boxscore_players_distinct$url[i])$shootsCatches
}

boxscore_players_distinct <- boxscore_players_distinct %>% select(c(playerId, fullname, shootsCatches, team))

faceoffs <- faceoffs %>% left_join((boxscore_players_distinct %>%
                                      mutate(playerId = as.numeric(playerId)) %>%
                                     rename("p1_fullname" = "fullname", "p1_team" = "team",
                                            "p1_shoots" = shootsCatches)),
                                   by = c("player_1" = "playerId"))

faceoffs <- faceoffs %>% left_join((boxscore_players_distinct %>%
                                      mutate(playerId = as.numeric(playerId)) %>%
                                      rename("p2_fullname" = "fullname", "p2_team" = "team",
                                             "p2_shoots" = shootsCatches)),
                                   by = c("player_2" = "playerId"))

faceoffs <- faceoffs %>% mutate(p1_faceoff_side = yCoord/22 * ifelse(p1_team == homeTeam, 1, -1) *
                                  ifelse(homeTeamDefendingSide == "right", 1, -1) * 
                                  ifelse(p1_shoots == "L", -1, 1) + 2) %>% 
  mutate(p2_faceoff_side = yCoord/22 * ifelse(p2_team == homeTeam, 1, -1) * 
           ifelse(homeTeamDefendingSide == "right", 1, -1) * ifelse(p2_shoots == "L", -1, 1) + 2)


unique_faceoff_player <- faceoffs %>% select(c(player_1, player_2)) %>% 
  pivot_longer(cols = everything()) %>% group_by(value) %>% 
  summarize(n_wins = sum(name == "player_1"), n = n(), prop = n_wins/n) %>%
  arrange(desc(n)) %>% mutate(index = as.numeric(rownames(.)))

unique_faceoff_player <- unique_faceoff_player %>% 
  left_join((boxscore_players_distinct %>%
               mutate(playerId = as.numeric(playerId))), by = c("value" = "playerId"))

faceoff_matrix_maker <- faceoffs %>% select(c(player_1, player_2)) %>%
  left_join((unique_faceoff_player %>% select(c(value, index))), 
            by = c("player_1" = "value")) %>% rename("player_1_index" = "index") %>%
  left_join((unique_faceoff_player %>% select(c(value, index))), 
            by = c("player_2" = "value")) %>% rename("player_2_index" = "index") %>%
  select(c(player_1_index, player_2_index))

#The Bradley-Terry model can run into problems when entries in a head-to-head matchup 
#are not part of the connected network of all other entries. This could happen in our
#case if two players took a faceoff against each other and never against anyone else.
#While the nature of faceoffs makes this extremely improbable over a whole season, we
#verify just in case.

faceoff_matrix <- matrix(0, nrow = nrow(unique_faceoff_player) + 1, ncol = nrow(unique_faceoff_player) + 1)  

for (i in 1:nrow(faceoff_matrix_maker)){
  faceoff_matrix[faceoff_matrix_maker$player_1_index[i],faceoff_matrix_maker$player_2_index[i]] <- 
    faceoff_matrix[faceoff_matrix_maker$player_1_index[i],faceoff_matrix_maker$player_2_index[i]] + 1
}
for (i in 2:nrow(faceoff_matrix)){
  for (j in 1:(i-1)){
    faceoff_matrix[i,j] <- faceoff_matrix[i,j] + faceoff_matrix[j,i]
  }
}

unaccounted_indices <- 2:nrow(unique_faceoff_player)
connected_indices <- 1
added_indices <- list()
removed_rows <- 1
while (removed_rows > 0){
  for (i in unaccounted_indices){
    for (j in connected_indices){
      if ((faceoff_matrix[i,j]) > 0) {
        added_indices <- append(added_indices, i) %>% unique()
        unaccounted_indices <- setdiff(unaccounted_indices, i)
      }
    }
  }
  removed_rows <- length(added_indices) - length(connected_indices)
  connected_indices <- added_indices
}
#all accounted for

faceoff_data <- list("N" = nrow(faceoffs), "K" = nrow(unique_faceoff_player),
                     "y" = rep(as.integer(1),nrow(faceoffs)), 
                     "player_1" = as.integer(faceoff_matrix_maker$player_1_index),
                     "player_2" = as.integer(faceoff_matrix_maker$player_2_index))

#The first fit will just include player variables, not accounting 
#for handedness. We fit the model and simulate new observations
#(y_pred) from the parameters to test our method.

basic_bradley_terry_fit <- stan(model_code = "
data{
  int N;
  int K;
  array[N] int <lower=0, upper=1> y;
  array[N] int <lower=1, upper=K> player_1;
  array[N] int <lower=1, upper=K> player_2;
}                                
parameters{
  vector[K] theta;
}
model{
  theta ~ normal(0, 0.33);
  y ~ bernoulli_logit(theta[player_1] - theta[player_2]);
}
generated quantities{
  array[N] int y_pred = bernoulli_logit_rng(theta[player_1] - theta[player_2]);
}
", data = faceoff_data, cores = mc.cores)

faceoff_samples <- rstan::extract(basic_bradley_terry_fit)

faceoff_quantiles <- apply(faceoff_samples$theta, 2, function(x) quantile(x, probs = seq(0.1,0.9,0.1)))

unique_faceoff_player %>% bind_cols(t(faceoff_quantiles)) %>% view()
#the coefficients can be interpreted as the proportional change in log odds
#of succcess when a given player is taking part in a faceoff. Negative coefficients
#are bad and positive coefficients are good. The Bayesian approach allows us to 
#easily present the uncertainty in parameter samples as quantiles. The 50% quantile
#for instance represents the median parameter sample. This is valuable in a case 
#like this where players have vastly different numbers of faceoffs taken.

#We could easily represent these coefficients as the percentage of faceoffs won
#given base odds of probability of success/probability of failure = 1/1 (i.e.
#50% chance of success):

a <- exp(faceoff_samples$theta)/(1 +  exp(faceoff_samples$theta)) 
faceoff_percent_quantiles <- apply(a, 2, function(x) quantile(x, probs = seq(0.1,0.9,0.1)))

#we see the necessity of modeling if we compare the difference between
#the empirical percentage of faceoffs won and the percentage won implied
#by our model:

unique_faceoff_percents <- unique_faceoff_player %>% bind_cols(t(faceoff_percent_quantiles)) %>%
  mutate(diff = `50%` - prop)

#observing the diff column here shows dramatic differences once we account
#for just the number of faceoffs taken leading to far greater uncertainty
#and potential large differences in opponent quality.

plot(unique_faceoff_percents$n, unique_faceoff_percents$diff, pch = 20, cex = 0.1)

#Not only does this plot show as expected that players with fewer than 50 faceoffs
#taken are exremely noisy, but there's a clear but very slight positive slope in
#difference for players with increasing faceoffs taken. Players who are good at faceoffs
#take a lot more of them, and face each other more often, so that the percentage won 
#underrates the true strength of these players against the average player.


#we can consider various posterior retrodictive summaries to compare how well
#our observations simulated from theta parameter samples compare to the observed
#data.
y_post_pred <- faceoff_samples$y_pred
sum(y_post_pred)/(4000 * nrow(faceoffs)) #.512
player_wins <- bind_cols((faceoffs %>% select(c(player_1, player_2))),
          "player_1_wins" = apply(y_post_pred, 2, 
                                  function(x) sum(x))) %>%
  mutate("player_2_wins" = 4000 - player_1_wins)
player_wins <- player_wins %>% pivot_longer(cols = c(player_1_wins, player_2_wins)) %>%
  mutate(player = ifelse(name == "player_1_wins", player_1, player_2)) %>%
  select(-c(player_1, player_2))

player_wins_summary <- player_wins %>% group_by(player) %>%
  summarize(n_pred_wins = sum(value)/4000) %>% left_join(unique_faceoff_player %>%
                                            select(c(fullname, value, n_wins, n, prop)), by =
                                            c("player" = "value"))
#the raw residuals seem very reasonable
hist((player_wins_summary$n_wins - player_wins_summary$n_pred_wins), breaks = 30)

#however we can see a linear pattern in the residuals: players
#who take a lot of faceoffs are slightly but consistently
#underestimated. 

plot(player_wins_summary$n, (player_wins_summary$n_wins - player_wins_summary$n_pred_wins), 
     pch = 20, cex = 0.1)

#I suspect that this is due to the fact that our prior
#assumes a neutral value for a player's strength at faceoffs,
#but the data is not informative for players who take very few 
#faceoffs, who are also presumably weaker at them, or the coach would
#trust them to take more. Inferences for these players cannot be pulled
#away from the neutral prior by the weak data, so the model overestimates
#them, and when they go against a strong player, our model underestimates
#that strong player's chances. Even though the proportion of residual 
#error is small for strong players with lots of data points (maybe 5
#predicted wins fewer out of 1500+ faceoffs), we should still try to 
#adjust our model.

#We can start by simply putting a lower prior on theta. The 0.33 value
#for the standard deviation implies a success rate of about 8%. If we
#simply make the mean of our theta variable -0.33, implying a success
#rate of around 42%, we can hope that our data will resolve players with
#a lot of data close to their true value while not overestimating players
#who only take a few faceoffs.

bradley_terry_fit_lower_prior <- stan(model_code = "
data{
  int N;
  int K;
  array[N] int <lower=0, upper=1> y;
  array[N] int <lower=1, upper=K> player_1;
  array[N] int <lower=1, upper=K> player_2;
}                                
parameters{
  vector[K] theta;
}
model{
  theta ~ normal(-0.33, 0.33);
  y ~ bernoulli_logit(theta[player_1] - theta[player_2]);
}
generated quantities{
  array[N] int y_pred = bernoulli_logit_rng(theta[player_1] - theta[player_2]);
}
", data = faceoff_data, cores = mc.cores)

faceoff_samples_lower_prior <- rstan::extract(bradley_terry_fit_lower_prior)

#we can consider various posterior retrodictive summaries to compare how well
#our observations simulated from theta parameter samples compare to the observed
#data.
y_post_pred_lower_prior <- faceoff_samples_lower_prior$y_pred
sum(y_post_pred_lower_prior)/(4000 * nrow(faceoffs)) #.512
player_wins_lower_prior <- bind_cols((faceoffs %>% select(c(player_1, player_2))),
                         "player_1_wins" = apply(y_post_pred_lower_prior, 2, 
                                                 function(x) sum(x))) %>%
  mutate("player_2_wins" = 4000 - player_1_wins)
player_wins_lower_prior <- player_wins_lower_prior %>%
  pivot_longer(cols = c(player_1_wins, player_2_wins)) %>%
  mutate(player = ifelse(name == "player_1_wins", player_1, player_2)) %>%
  select(-c(player_1, player_2))

player_wins_summary_lower_prior <- player_wins_lower_prior %>%
  group_by(player) %>% summarize(n_pred_wins = sum(value)/4000) %>%
  left_join(unique_faceoff_player %>% 
              select(c(fullname, value, n_wins, n, prop)), by =
              c("player" = "value"))


hist((player_wins_summary_lower_prior$n_wins - player_wins_summary_lower_prior$n_pred_wins), breaks = 30)

#the pattern persists for the most part

plot(player_wins_summary_lower_prior$n, (player_wins_summary_lower_prior$n_wins - player_wins_summary_lower_prior$n_pred_wins), 
     pch = 20, cex = 0.1)

#we alter our model so that our prior for our player coefficient is conditional
#on the number of faceoffs taken. Considering that we see a linear pattern in 
#the residuals, a normal model for n_faceoffs taken may be appropriate:

player_wins_summary <- player_wins_summary %>% arrange(desc(n)) %>%
  mutate(scaled_n = (n - (max(n)/2))/(max(n)/2))

faceoff_data_conditional <- list(
  "N" = nrow(faceoffs), "K" = nrow(unique_faceoff_player),
  "y" = rep(as.integer(1),nrow(faceoffs)), 
  "player_1" = as.integer(faceoff_matrix_maker$player_1_index),
  "player_2" = as.integer(faceoff_matrix_maker$player_2_index),
  "N_taken" = player_wins_summary$scaled_n
)

bradley_terry_sim_conditional <- stan(model_code = "
data{
  int N;
  int K;
  array[N] int <lower=0, upper=1> y;
  array[N] int <lower=1, upper=K> player_1;
  array[N] int <lower=1, upper=K> player_2;
  vector<lower=-1, upper=1>[K] N_taken;
}                                
transformed data{
  real beta = normal_rng(0.33,0.1);
  vector[K] theta_prior = beta*N_taken;
}
generated quantities{
  vector[K] theta;
  for (k in 1:K){
    theta[k] = normal_rng(theta_prior[k], 0.1); 
  }
  array[N] int y_pred = bernoulli_logit_rng(theta[player_1] - theta[player_2]);
}
", data = faceoff_data_conditional, cores = mc.cores, chains = 1, iter = 1, warmup = 0,
algorithm = "Fixed_param")

conditional_sim <- rstan::extract(bradley_terry_sim_conditional)
plot(player_wins_summary$n, conditional_sim$theta, pch = 20, cex = 0.1)

bradley_terry_fit_conditional <- stan(model_code = "
data{
  int N;
  int K;
  array[N] int <lower=0, upper=1> y;
  array[N] int <lower=1, upper=K> player_1;
  array[N] int <lower=1, upper=K> player_2;
  vector<lower=-1, upper=1>[K] N_taken;
}                                
parameters{
  vector[K] theta;
  real beta;
}
transformed parameters{
  vector[K] theta_prior = beta*N_taken;
}
model{
  beta ~ normal(0.33,0.1); 
  theta ~ normal(theta_prior, 0.33); 
  y ~ bernoulli_logit(theta[player_1] - theta[player_2]);
}
generated quantities{
  array[N] int y_pred = bernoulli_logit_rng(theta[player_1] - theta[player_2]);
}
", data = faceoff_data_conditional, cores = mc.cores)

faceoff_samples_conditional <- rstan::extract(bradley_terry_fit_conditional)
rm(bradley_terry_fit_conditional)

y_post_pred_conditional <- faceoff_samples_conditional$y_pred
sum(y_post_pred_conditional)/(4000 * nrow(faceoffs)) #.513
player_wins_conditional <- bind_cols((faceoffs %>% select(c(player_1, player_2))),
                                     "player_1_wins" = apply(y_post_pred_conditional, 2, 
                                                             function(x) sum(x))) %>%
  mutate("player_2_wins" = 4000 - player_1_wins)
player_wins_conditional <- player_wins_conditional %>%
  pivot_longer(cols = c(player_1_wins, player_2_wins)) %>%
  mutate(player = ifelse(name == "player_1_wins", player_1, player_2)) %>%
  select(-c(player_1, player_2))

player_wins_summary_conditional <- player_wins_conditional %>%
  group_by(player) %>% summarize(n_pred_wins = sum(value)/4000) %>%
  left_join(unique_faceoff_player %>% 
              select(c(fullname, value, n_wins, n, prop)), by =
              c("player" = "value"))


hist((player_wins_summary_conditional$n_wins - player_wins_summary_conditional$n_pred_wins), breaks = 30)

#the residuals are in the same range but no longer have a pattern of underestimating
#the best players:

plot(player_wins_summary_conditional$n, (player_wins_summary_conditional$n_wins - player_wins_summary_conditional$n_pred_wins), 
     pch = 20, cex = 0.1)

#the parameters are no longer centered at zero due to our prior
cond_mean <- mean(faceoff_samples_conditional$theta)
cond_mean
#[1] -0.2811865

hist(apply(faceoff_samples$theta, 2, function(x) mean(x)))
hist(apply(faceoff_samples_conditional$theta, 2, function(x) mean(x)))

#recenter the parameters around the mean strength
faceoff_quantiles_conditional <- apply(faceoff_samples_conditional$theta - cond_mean, 2, function(x) quantile(x, probs = seq(0.1,0.9,0.1)))

unique_faceoff_player %>% bind_cols(t(faceoff_quantiles_conditional)) %>% view()

a3 <- exp(faceoff_samples_conditional$theta - cond_mean)/(1 +  exp(faceoff_samples_conditional$theta - cond_mean)) 
faceoff_percent_quantiles_conditional <- apply(a3, 2, function(x) quantile(x, probs = seq(0.1,0.9,0.1)))

#we see the necessity of modeling if we compare the difference between
#the empirical percentage of faceoffs won and the percentage won implied
#by our model:

unique_faceoff_percents_conditional <- unique_faceoff_player %>% bind_cols(t(faceoff_percent_quantiles_conditional)) %>%
  mutate(diff = `50%` - prop)

#under this model, the strongest players appear even stronger,
#and have a better chance of winning a faceoff against the perfect
#average player than under the previous model that underrated them.

unique_faceoff_conditional_percents <- unique_faceoff_player %>%
  bind_cols(t(faceoff_percent_quantiles[c(1,5,9),])) %>%
  rename_with(.cols = c(9:11), ~ paste(., "basic")) %>%
  bind_cols(t(faceoff_percent_quantiles_conditional[c(1,5,9),])) %>%
  rename_with(.cols = c(12:14), ~ paste(., "conditional")) %>% 
  mutate(diff2 = `50% conditional` - `50% basic`)

#plot the difference in implied probability of winning a draw against
#an average player from the basic model to the conditional prior
#model against the number of faceoffs taken 

plot(unique_faceoff_conditional_percents$n, 
     unique_faceoff_conditional_percents$diff2, pch = 20, cex = 0.1)

#Now we want to include whether a draw is on a player's strong or weak side to see
#how this affects the coefficients. We include a parameter for the number of 
#faceoffs a player has taken and make this the prior for theta as in the last example.


faceoff_data_side <- list("N" = nrow(faceoffs), "K" = nrow(unique_faceoff_player),
                     "y" = rep(as.integer(1),nrow(faceoffs)), 
                     "player_1" = as.integer(faceoff_matrix_maker$player_1_index),
                     "player_2" = as.integer(faceoff_matrix_maker$player_2_index),
                     "player_1_side" = as.integer(faceoffs$p1_faceoff_side),
                     "player_2_side" = as.integer(faceoffs$p2_faceoff_side),
                     "N_taken" = faceoff_data_conditional$N_taken)

bradley_terry_side_fit <- stan(model_code = "
data{
  int N;
  int K;
  array[N] int <lower=0, upper=1> y;
  array[N] int <lower=1, upper=K> player_1;
  array[N] int <lower=1, upper=K> player_2;
  array[N] int <lower=1, upper=3> player_1_side;
  array[N] int <lower=1, upper=3> player_2_side;
  vector<lower=-1, upper=1>[K] N_taken;
}                                
parameters{
  vector[K] theta;
  vector[2] alpha;
  real beta;
}
transformed parameters{
  vector[3] alpha_trans = [alpha[1], 0, alpha[2]]';
  vector[K] theta_prior = beta*N_taken;
}
model{
  beta ~ normal(0.33,0.1); 
  theta ~ normal(theta_prior, 0.33); 
  alpha ~ normal(0, 0.33);
  y ~ bernoulli_logit((theta[player_1] + alpha_trans[player_1_side]) - (theta[player_2] + alpha_trans[player_2_side]));
}
generated quantities{
  array[N] int y_pred = bernoulli_logit_rng((theta[player_1] + alpha_trans[player_1_side]) - (theta[player_2] + alpha_trans[player_2_side]));
}
", data = faceoff_data_side, cores = mc.cores)

bradley_terry_side_samples <- rstan::extract(bradley_terry_side_fit)
rm(bradley_terry_side_fit)
view(apply(bradley_terry_side_samples$alpha, 2, function(x) quantile(x, probs = seq(0.1,0.9,0.1))))
#the quantiles are very wide but median on side value gives you a change in odds of success of
#exp(.13) = 1.16, equivalent to 53.2% probability of winning an otherwise 50-50 faceoff.

#as before, the samples are not centered anymore
faceoff_quantiles_side <- apply(bradley_terry_side_samples$theta, 2, function(x) quantile(x, probs = seq(0.1,0.9,0.1)))
side_mean <- mean(bradley_terry_side_samples$theta)

a2 <- exp(bradley_terry_side_samples$theta - side_mean)/(1 +  exp(bradley_terry_side_samples$theta - side_mean)) 
faceoff_percent_quantiles_side <- apply(a2, 2, function(x) quantile(x, probs = seq(0.1,0.9,0.1)))

unique_faceoff_side_percents <- unique_faceoff_player %>%
  bind_cols(t(faceoff_percent_quantiles_conditional[c(1,5,9),])) %>%
  rename_with(.cols = c(9:11), ~ paste(., "no_side")) %>%
  bind_cols(t(faceoff_percent_quantiles_side[c(1,5,9),])) %>%
  rename_with(.cols = c(12:14), ~ paste(., "with_side")) %>% 
  mutate(diff2 = `50% with_side` - `50% no_side`)

mean(abs(unique_faceoff_side_percents$diff2))

hist(unique_faceoff_side_percents$diff2, breaks = 30)
#most players only see a very slight change in their faceoff strength
#with only a few seeing a change in over 1%

#no strong pattern here:
plot(unique_faceoff_side_percents$n, unique_faceoff_side_percents$diff2,
     pch = 20, cex = 0.1)

y_post_pred_side <- bradley_terry_side_samples$y_pred

#the model is making slightly more correct predictions
#than the previous model  
sum(y_post_pred_side)/(4000 * nrow(faceoffs)) #.516
player_wins_side <- bind_cols((faceoffs %>% select(c(player_1, player_2))),
                                     "player_1_wins" = apply(y_post_pred_side, 2, 
                                                             function(x) sum(x))) %>%
  mutate("player_2_wins" = 4000 - player_1_wins)
player_wins_side <- player_wins_side %>%
  pivot_longer(cols = c(player_1_wins, player_2_wins)) %>%
  mutate(player = ifelse(name == "player_1_wins", player_1, player_2)) %>%
  select(-c(player_1, player_2))

player_wins_summary_side <- player_wins_side %>%
  group_by(player) %>% summarize(n_pred_wins = sum(value)/4000) %>%
  left_join(unique_faceoff_player %>% 
              select(c(fullname, value, n_wins, n, prop)), by =
              c("player" = "value"))

#looks quite similar to conditional model plot
hist((player_wins_summary_side$n_wins - player_wins_summary_side$n_pred_wins), breaks = 30)

#again quite similar
plot(player_wins_summary_side$n, (player_wins_summary_side$n_wins - player_wins_summary_side$n_pred_wins), 
     pch = 20, cex = 0.1)

Raw_vs_modeled_comparison <-  
  bind_cols((unique_faceoff_player %>% select(c(fullname, prop))), 
          t(faceoff_percent_quantiles_side)) %>% arrange(`50%`)

Raw_vs_modeled_comparison$index <- as.numeric(1:nrow(Raw_vs_modeled_comparison))
Raw_vs_modeled_comparison <-   as.data.frame(apply(Raw_vs_modeled_comparison, 2,
                                     function(x) rep(x, each = 2)))
Raw_vs_modeled_comparison[,c(2:12)] <- 
  apply(Raw_vs_modeled_comparison[,c(2:12)], 2, function(x) as.numeric(x))

Raw_vs_modeled_comparison <- 
  Raw_vs_modeled_comparison %>% mutate(
    index = ifelse(as.numeric(rownames(.)) %% 2 == 0, index + 0.5, index - 0.5))

par(mar=c(4, 1, 1, 3))

plot(1, type = "n", xlim = c(0.4, 0.8), ylim = c(518,566),
      main = "Top Faceoff Players in 23/24 Season", 
     xlab = "Implied Neutral Win Percentage", yaxt = "n")
polygon(c(Raw_vs_modeled_comparison$`10%`,rev(Raw_vs_modeled_comparison$`90%`)),
      c(Raw_vs_modeled_comparison$index,rev(Raw_vs_modeled_comparison$index)),
      col = alpha("coral", 0.2), border = NA)
polygon(c(Raw_vs_modeled_comparison$`20%`,rev(Raw_vs_modeled_comparison$`80%`)),
        c(Raw_vs_modeled_comparison$index,rev(Raw_vs_modeled_comparison$index)),
        col = alpha("coral", 0.4), border = NA)
polygon(c(Raw_vs_modeled_comparison$`30%`,rev(Raw_vs_modeled_comparison$`70%`)),
        c(Raw_vs_modeled_comparison$index,rev(Raw_vs_modeled_comparison$index)),
        col = alpha("coral", 0.6), border = NA)
polygon(c(Raw_vs_modeled_comparison$`40%`,rev(Raw_vs_modeled_comparison$`60%`)),
        c(Raw_vs_modeled_comparison$index,rev(Raw_vs_modeled_comparison$index)),
        col = alpha("coral", 0.8), border = NA)
lines(Raw_vs_modeled_comparison$`50%`, Raw_vs_modeled_comparison$index,
      col = "red")

plot_labels <- 
  Raw_vs_modeled_comparison$fullname[as.numeric(
    rownames(Raw_vs_modeled_comparison)) %% 2 == 0]
for (i in 517:566){
  text(x = 0.5, y = i, label = plot_labels[i], cex = 0.5)
}


