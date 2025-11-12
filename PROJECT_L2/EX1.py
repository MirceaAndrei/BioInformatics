from itertools import product

def generate_combinations(alphabet, length):
    """Generate all combinations of given alphabet with specified length."""
    return [''.join(p) for p in product(alphabet, repeat=length)]


def count_occurrences(text, pattern):
    """Count occurrences of pattern in text (case-insensitive)."""
    text = text.upper()
    pattern = pattern.upper()
    count = 0
    n, m = len(text), len(pattern)
    
    for i in range(n - m + 1):
        if text[i:i+m] == pattern:
            count += 1
    return count


def main():
    input_str = input("Enter a string: ").strip()

    # Extract alphabet (unique letters)
    alphabet = sorted(set(c.upper() for c in input_str if c.isalpha()))
    print("Alphabet of string:", ''.join(alphabet))

    # Generate 2- and 3-character combinations
    combos2 = generate_combinations(alphabet, 2)
    combos3 = generate_combinations(alphabet, 3)

    total_pairs = max(0, len(input_str) - 1)
    print("\n--- 2-character combinations ---")
    for combo in combos2:
        count = count_occurrences(input_str, combo)
        percentage = (100.0 * count / total_pairs) if total_pairs > 0 else 0.0
        print(f"{combo}: {percentage:.2f}% ({count} occurrences)")

    total_triplets = max(0, len(input_str) - 2)
    print("\n--- 3-character combinations ---")
    for combo in combos3:
        count = count_occurrences(input_str, combo)
        percentage = (100.0 * count / total_triplets) if total_triplets > 0 else 0.0
        print(f"{combo}: {percentage:.2f}% ({count} occurrences)")


if __name__ == "__main__":
    main()
