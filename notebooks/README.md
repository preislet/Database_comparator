# 📝 Prepared Notebooks for the Project

If you're using the Docker package, the notebooks provided in the `notebooks` folder are ready to use without any modifications.

However, if you're running the code from source and haven't installed the database-comparator package, you’ll need to manually adjust the import in the notebooks. Specifically, replace the first cell:

```python
import pandas as pd
import numpy as np

from Database_comparator import db_compare
```

with the following workaround:

```python
import pandas as pd
import numpy as np
import os
import sys
# Get the absolute path of the project root directory
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# open folder Database_comparator
sys.path.append(project_root)

# Now you can import db_compare
import Database_comparator.db_compare as db_compare
```

## 🔧 Important Note

This modification will work only if youac run the notebooks from the `notebooks` folder. If you run them from another location, you will need to adjust the path accordingly or install the package properly.

## ⚠️ Recommendation

To avoid setup issues and ensure a smoother experience, we strongly recommend using the Docker image provided with this project.