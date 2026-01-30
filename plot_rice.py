import pandas as pd
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
# 1. Names of your 3 CSV files
base_path = "csv/"
file_paths = [f"{base_path}japonais_cleaned.csv", f"{base_path}camargue_cleaned.csv", f"{base_path}basmati_cleaned.csv"]

# 2. Names to appear on the X-axis
labels = ["Japonais", "Camargue", "Basmati"]

# 3. Features you want to compare (Must match CSV headers)
features_to_plot = [
    {"col": "PolyArea",      "title": "Polygon Area Distribution",      "ylabel": "Area (pixels)"},
    {"col": "PolyPerimeter", "title": "Polygon Perimeter Distribution", "ylabel": "Perimeter (pixels)"},
    {"col": "Circularity",   "title": "Circularity Distribution",       "ylabel": "Circularity (0-1)"}
]

def main():
    # 1. Load Data
    data_frames = []
    for f in file_paths:
        try:
            df = pd.read_csv(f)
            data_frames.append(df)
            print(f"Loaded {f}: {len(df)} grains.")
        except FileNotFoundError:
            print(f"Error: Could not find file '{f}'. Make sure it is in the same folder.")
            return

    # 2. Setup Plot (1 Row, 3 Columns)
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))

    # 3. Plot Each Feature
    for i, feature in enumerate(features_to_plot):
        ax = axes[i]
        col_name = feature["col"]
        
        # Collect data for this specific column from all 3 files
        # We use .dropna() to be safe against empty rows
        data_to_plot = [df[col_name].dropna() for df in data_frames]
        
        # Create Boxplot
        # patch_artist=True allows us to fill them with color
        bplot = ax.boxplot(data_to_plot, 
                           labels=labels,
                           patch_artist=True,
                           medianprops=dict(color="black", linewidth=1.5))
        
        # Formatting
        ax.set_title(feature["title"], fontsize=14, fontweight='bold')
        ax.set_ylabel(feature["ylabel"], fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.6)

        # Optional: Color the boxes for visual distinction
        colors = ['lightblue', 'lightgreen', 'lightpink']
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)

    # 4. Final Layout Adjustments
    plt.tight_layout()
    
    # 5. Save and Show
    output_filename = "rice_comparison_boxplots.png"
    plt.savefig(output_filename, dpi=300)
    print(f"Plot saved to {output_filename}")
    plt.show()

if __name__ == "__main__":
    main()