import csv

if __name__ == '__main__':
    with open("test.csv", encoding='utf-8') as r_file:
        data = csv.DictReader(r_file, delimiter=',')
        for row in data:
            print(row["x"], row['f'])
