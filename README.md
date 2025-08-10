# data_and_info

# Dataset Catalog Specification

This document describes the expected structure of the **dataset catalog JSON** used by extracting data in API systems.

The catalog defines curated datasets, their metadata, and downloadable artifacts for local use.

---

## File Format

* The catalog is a **single JSON object** with the following **top-level keys**:

| Key        | Type     | Required | Description                                                                                                       |
| ---------- | -------- | -------- | ----------------------------------------------------------------------------------------------------------------- |
| `base_url` | `string` | No       | Optional default base URL for resolving `artifact.path` when downloading datasets. Can be overridden in `.get()`. |
| `datasets` | `array`  | Yes      | List of dataset entries (see [Dataset Object](#dataset-object) below).                                            |

---

## Dataset Object

Each entry in `datasets` is a JSON object with the following fields:

| Field         | Type     | Required | Description                                                                                                |
| ------------- | -------- | -------- | ---------------------------------------------------------------------------------------------------------- |
| `name`        | `string` | Yes      | Unique dataset identifier (short, lowercase, underscores).                                                 |
| `title`       | `string` | No       | Human-readable dataset name. Defaults to `name` if missing.                                                |
| `version`     | `string` | Yes      | Dataset version (e.g., `"1.0.0"` or `"2025-08"`).                                                          |
| `license`     | `string` | No       | License string (e.g., `"CC-BY-4.0"`, `"MIT"`, `"Proprietary"`). `"NA"` if unspecified.                     |
| `description` | `string` | No       | Short description of the dataset content and purpose.                                                      |
| `tasks`       | `array`  | No       | List of tasks/use cases this dataset supports (e.g., `["enzyme_classification", "structure_clustering"]`). |
| `artifacts`   | `array`  | Yes      | List of downloadable files (see [Artifact Object](#artifact-object) below).                                |
| `homepage`    | `string` | No       | URL to dataset homepage or documentation.                                                                  |
| `extras`      | `object` | No       | Arbitrary additional metadata (key/value pairs).                                                           |

---

## Artifact Object

Each entry in the `artifacts` array is a JSON object with:

| Field      | Type      | Required | Description                                                                                                            |
| ---------- | --------- | -------- | ---------------------------------------------------------------------------------------------------------------------- |
| `path`     | `string`  | Yes      | Relative path of the artifact to be appended to `base_url`.                                                            |
| `sha256`   | `string`  | Yes      | SHA-256 checksum of the file (hex string). Used for integrity verification.                                            |
| `size`     | `integer` | Yes      | File size in bytes. Used for progress bars and validation.                                                             |
| `optional` | `boolean` | No       | If `true`, this artifact is downloaded only if explicitly requested with `include_optional=True`. Defaults to `false`. |

---

## Example

```json
{
  "base_url": "https://example.org/datasets",
  "datasets": [
    {
      "name": "enzymes_demo",
      "title": "Demo Enzyme Classification Dataset",
      "version": "1.0.0",
      "license": "CC-BY-4.0",
      "description": "A small set of protein sequences and annotations for testing enzyme classification workflows.",
      "tasks": ["enzyme_classification", "benchmark"],
      "artifacts": [
        {
          "path": "enzymes_demo_v1.csv",
          "sha256": "5e884898da28047151d0e56f8dc6292773603d0d6aabbddc1b8e9f4e6b82e7aa",
          "size": 123456,
          "optional": false
        },
        {
          "path": "enzymes_demo_metadata.json",
          "sha256": "f4b645a1b02ff17daaa2c7f3f019bcfc0ec9b02eeaf8ed22ad22df7a1f48d0f2",
          "size": 2345,
          "optional": true
        }
      ],
      "homepage": "https://example.org/enzymes_demo",
      "extras": {
        "num_sequences": 1000,
        "source": "UniProtKB",
        "date_created": "2025-08-01"
      }
    }
  ]
}
```

---

## Notes

1. **Base URL resolution**

   * If `base_url` exists at the top level, it will be **prepended** to all artifact `path`s unless overridden by `DatasetRegistry.get(base_url=...)`.
   * If no `base_url` is present, you **must** provide one when calling `.get()`.

2. **Multiple versions**

   * You may have multiple entries with the same `name` but different `version` values.
   * `DatasetRegistry.info(name)` returns the **latest** version (last in the list) if no version is specified.

3. **Optional artifacts**

   * Marking `optional: true` allows users to skip large or auxiliary files unless they request them.

4. **Checksum verification**

   * The SHA-256 hash is verified after download to ensure file integrity.

