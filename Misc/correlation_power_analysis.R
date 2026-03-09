library(dplyr)
library(ggplot2)

prediction_test <- function(df, n_sample, n_iter) {
  RMSD = rep(0, times= n_iter)
  for (i_var in seq(from=1, to=n_iter)) {
    test_group = df[sample(seq(from=1, to=n_common), n_sample, replace = FALSE),]
    model_fit = lm(set2 ~ poly(set1, degree = 4), data = test_group)
    projected_value = predict(model_fit, df)
    RMSD[i_var] = sqrt(mean((projected_value - df$set2)^2))
  }
  return(mean(RMSD))
}

data_set1 = read.table("")
data_set1 = rename(data_set1, set1 = x)
data_set1 = mutate(data_set1, set1 = 2^set1)
data_set2 = read.table("")
data_set2 = rename(data_set2, set2 = x)
data_set2 = mutate(data_set2, set2 = 2^set2)

common_dataset = left_join(data_set1, data_set2, by = "Lysate.code")
common_dataset = filter(common_dataset, !(is.na(set1)|is.na(set2)))

n_common = nrow(common_dataset)
min_RMSD = prediction_test(common_dataset, n_common, 1)
i_size = 10
threshold = 1.1 #multiplier similarity approached from the top
score = c()
size = c()

while (i_size <= n_common) {
  mean_score = prediction_test(common_dataset, i_size, 5)
  score = append(score, mean_score)
  size = append(size, i_size)
  
  if (mean_score <= min_RMSD * threshold | i_size == n_common) {
    local_ceiling = i_size
    i_size = prev_i
    break()
  } else {
    prev_i = i_size
    i_size = min(i_size*2, n_common) 
  }
}

while (i_size < local_ceiling) {
  mid = (i_size + local_ceiling) %/% 2
  mean_score = prediction_test(common_dataset, mid, 5)
  score = append(score, mean_score)
  size = append(size, mid)
  
  if (mean_score <= min_RMSD * threshold) {
    local_ceiling = mid
  } else {
    i_size = mid +1 
  }
}

print(paste((threshold-1)*100, "% threshold reached at ", mid, "/", n_common, sep = ""))
#plot(size[c(-1, -2, -3, -4)], score[c(-1, -2, -3, -4)], type = "b")
curve = data.frame( set1 = seq(from = min(common_dataset$set1), to = max(common_dataset$set1), length.out =20))
curvem = lm(set2 ~ poly(set1, degree = 4), data = common_dataset[sample(seq(from=1, to=n_common), mid, replace = FALSE),])
curve$set2 = predict(curvem, curve)
corr_plot = ggplot() +
  geom_point(data = common_dataset, aes(x = set1, y = set2)) +
  geom_line(data = curve, aes( x = set1, y = set2)) +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2")
plot(corr_plot)