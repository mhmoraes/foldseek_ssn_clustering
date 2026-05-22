# ==============================================================================
# Custom Graphics Tweaking
# ==============================================================================
library(tidyverse)
library(ggraph)
library(tidygraph)
library(colorspace)

# 1. Instantly load the complete, saved network package back into memory
ssn_data <- readRDS("results/ssn_graph_object.rds")

network_graph     <- ssn_data$graph
network_coords    <- ssn_data$coordinates
optimized_cutoff  <- ssn_data$cutoff

# 2. Tweak your graphics rules cleanly!
# Example: Let's re-render the network map with custom parameters
tweaked_plot <- ggraph(network_graph, layout = "manual", x = network_coords[,1], y = network_coords[,2]) +
  # Make alignment lines subtly darker
  geom_edge_diagonal0(aes(alpha = lddt), color = "grey30", show.legend = FALSE) + 
  # Customize node points (Color nodes by community id, size nodes by node degree)
  geom_node_point(aes(color = community_id, size = node_degree), alpha = 0.90) +
  
  # CUSTOM TWEAK: Apply a striking dark palette or alternative colors pace profile
  scale_color_manual(values = colorspace::sequential_hcl(length(unique(V(network_graph)$community_id)), "Viridis")) +
  scale_size_continuous(name = "Node Degree (Hubness)", range = c(1.0, 8.0)) +
  
  theme_graph() +
  labs(title = "Custom Tweaked Structural Similarity Network",
       subtitle = sprintf("Interactive Modification Session | Fixed Cutoff lDDT >= %s", optimized_cutoff))

# Display the newly modified image instantly inside your RStudio Plots panel
print(tweaked_plot)