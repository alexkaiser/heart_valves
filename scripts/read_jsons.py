import json
import os

def find_positions(data):
    positions = []

    if isinstance(data, dict):
        for key, value in data.items():
            if key == "position":
                if isinstance(value, list):
                    positions.append(value)  # keep as a group
                else:
                    positions.append([value])  # wrap single value as list
            positions.extend(find_positions(value))
    elif isinstance(data, list):
        for item in data:
            positions.extend(find_positions(item))

    return positions


def main():
    for filename in os.listdir("."):
        if filename.startswith("aortic_annulus") and filename.endswith(".json"):
            file_path = os.path.join(".", filename)

            try:
                with open(file_path, "r", encoding="utf-8") as f:
                    data = json.load(f)

                positions = find_positions(data)
                name_without_ext = filename.split('.')[0] #os.path.splitext(filename)[0]

                print(f"{name_without_ext} = 0.1 * [")

                for group in positions:
                    line = ", ".join(str(p) for p in group)
                    print(line)

                print("]';\n")

            except Exception as e:
                print(f"Error processing {filename}: {e}\n")


if __name__ == "__main__":
    main()