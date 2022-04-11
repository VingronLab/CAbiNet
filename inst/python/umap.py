import umap

def python_umap(dm, n_neighbors, metric):
                  
  reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        metric=metric
        )
        
  embedding = reducer.fit_transform(dm)
  return embedding
