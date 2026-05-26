library(impower133)

dir.create("test_output", showWarnings = FALSE)

# 1. stimulate baseline data
cat("1/6 Simulating baseline data...\n")
df <- simulate_impower133(
  seed = 21,
  os_path = "overall_km/overall_IPD.csv",
  method = "random",
  output_path = "test_output/simulated_data.csv"
)

# 2. Table 1
cat("2/6 Generating Table 1...\n")
t1 <- make_table1(df)
gt::gtsave(t1, "test_output/table1.html")

# 3. Table 2
cat("3/6 Generating Table 2...\n")
t2 <- make_table2(df)
gt::gtsave(t2, "test_output/table2.html")

# 4. Table 3
cat("4/6 Generating Table 3...\n")
t3 <- make_table3(df)
gt::gtsave(t3, "test_output/table3.html")

# 5. OS K-M curve
cat("5/6 Plotting OS K-M curve...\n")
plot_survival("overall_km/overall_IPD.csv", type = "OS", output_path = "test_output/figure2a.png")

# 6. Forest plot
cat("6/6 Plotting forest plot...\n")
plot_figure2c(df, output_path = "test_output/figure2c.png")

cat("\nAll tests passed!\n")
