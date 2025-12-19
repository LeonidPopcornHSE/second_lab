import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("results.csv")
avg = df.groupby(["impl", "threads", "scenario"])["time"].mean().reset_index()

for scenario in avg["scenario"].unique():
    plt.figure()
    for impl in ["pthread", "my"]:
        sub = avg[(avg["scenario"] == scenario) & (avg["impl"] == impl)]
        plt.plot(sub["threads"], sub["time"], marker='o', label=impl)
    plt.title(f"Scenario: {scenario}")
    plt.xlabel("Threads")
    plt.ylabel("Time (seconds)")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"{scenario}.png")
    plt.close()
