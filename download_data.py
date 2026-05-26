import os
import crisprbrain

# 1. Initialize the official DataCommons pipeline client
client = crisprbrain.Client()

# 2. Assign the exact unique identifiers from your search screen
neuron_crispra_id = "Glutamatergic Neuron-RNA-Seq-CRISPRa-2020"
neuron_crispri_id = "Glutamatergic Neuron-RNA-Seq-CRISPRi-2020"
microglia_crispri_id = "iTF Microglia-Day-8-CROP-seq-CRISPRi"

# 3. Pull down the full structured screening matrix directly from the API cloud
print("Downloading global multi-omic perturbation matrices...")
neuron_activation_screen = client.screens[neuron_crispra_id]
neuron_inhibition_screen = client.screens[neuron_crispri_id]
microglia_inhibition_screen = client.screens[microglia_crispri_id]

# 4. Instantly convert the raw screen object into a Pandas DataFrame
df_neuron_a = neuron_activation_screen.to_data_frame()
df_neuron_i = neuron_inhibition_screen.to_data_frame()
df_microglia_i = microglia_inhibition_screen.to_data_frame()

# NEW STEP: Automatically create the folder if it doesn't exist yet!
output_dir = "./perturbation_data"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 5. Flush the matrices out into your local directory
df_neuron_a.to_csv(f"{output_dir}/neuron_crispra_DEG_matrix.tsv", sep="\t", index=False)
df_neuron_i.to_csv(f"{output_dir}/neuron_crispri_DEG_matrix.tsv", sep="\t", index=False)
df_microglia_i.to_csv(f"{output_dir}/microglia_crispri_DEG_matrix.tsv", sep="\t", index=False)

print("Matrices successfully exported to local pipeline workspace!")