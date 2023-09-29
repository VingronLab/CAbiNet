import umap.umap_ as umap

def python_umap(dm, n_neighbors, metric, seed=2358):
                  
  reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        metric=metric,
        random_state = seed
        )
        
  embedding = reducer.fit_transform(dm)
  return embedding
