#include "GA.hpp"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>

// インスタンス
// 各パラメータを初期化
GA::GA(const char* instance, int populationSize, int generationSize, float selectivePressure, float crossoverRatio, float crossoverProbability, float mutateProbability) : populationSize(populationSize), generationSize(generationSize), crossoverRatio(crossoverRatio), crossoverProbability(crossoverProbability), mutateProbability(mutateProbability)
{
    // rank-based法による選択確率の計算
    for (int i = populationSize; i > 0; --i)
    {
        float p = 2 - selectivePressure + 2 * (selectivePressure - 1) * (i - 1) / (populationSize - 1);
        selectProbabilities.push_back(p);
    }

    // インスタンスを外部ファイルからcitiesに読み込む
    ifstream ifs(instance);
    string line;

    while (getline(ifs, line))
    {
        istringstream iss(line);
        string iString, xString, yString;

        iss >> iString >> xString >> yString;

        float x = stof(xString);
        float y = stof(yString);

        City city(x, y);

        cities.push_back(city);
    }

    chromosomeSize = cities.size();

    srand((unsigned int)time(0));
}

// 初期集団の生成
// ランダムに順列を生成
void GA::initialize()
{
    best.fitness = 0;
    population.clear();

    for (int i = 0; i < populationSize; ++i)
    {
        Chromosome chromosome;
        population.push_back(chromosome);

        for (int j = 0; j < chromosomeSize; ++j)
        {
            population.at(i).genes.push_back(j);
        }

        int j = chromosomeSize;
        
        while (j > 0) {
            int k = rand() % j;
            j--;
            int t = population.at(i).genes.at(j);
            population.at(i).genes.at(j) = population.at(i).genes.at(k);
            population.at(i).genes.at(k) = t;
        }
    }
}

// 適合度の計算
// 経路の長さを適合度とする
void GA::evaluate()
{
    for (auto &individual : population)
    {
        individual.fitness = 0;
        for (int i = 0; i < individual.genes.size() - 1; ++i)
        {
            individual.fitness += dist(cities.at(individual.genes.at(i)), cities.at(individual.genes.at(i + 1)));
        }
        individual.fitness += dist(cities.at(individual.genes.at(individual.genes.size() - 1)), cities.at(individual.genes.at(0)));

        if (individual.fitness < best.fitness || best.fitness == 0)
        {
            best = individual;
        }
    }
}

// 選択
// 順位に対応した確率で選択を行う
void GA::select()
{
    matingPool.clear();

    sort(population.begin(), population.end());

    for (int i = 0; i < populationSize; ++i)
    {
        float x = (float)rand() / (float)RAND_MAX * (float)populationSize;
        float sum = 0.0f;

        for (int j = 0; j < populationSize; ++j)
        {
            sum += selectProbabilities.at(j);
            
            if (sum > x || j == populationSize - 1)
            {
                matingPool.push_back(population.at(j));
                break;
            }
        }
    }
}

// 交叉
// matingPoolのうちchromosomeRatioの割合だけ辺再組合せ交叉，残りを順序交叉で交叉を行う
void GA::crossover()
{
    int i = populationSize;

    while (i > 0) {
        int k = rand() % i;
        i--;
        Chromosome t = matingPool.at(i);
        matingPool.at(i) = matingPool.at(k);
        matingPool.at(k) = t;
    }

    population.clear();

    for (i = 0; i < populationSize / 2; ++i)
    {
        int parents[2] = { i * 2, i * 2 + 1 };

        float p = (float)rand() / (float)RAND_MAX;
        if (p < crossoverProbability)
        {
            // 辺再組合せ交叉
            if ((float)(i + 1) / (float)populationSize * 2.0f <= crossoverRatio)
            {
                vector<vector<int>> edgeTable(chromosomeSize, vector<int>());

                // 各都市に対する隣接リストの生成
                for (int j = 0; j < chromosomeSize; ++j)
                {
                    for (int k = 0; k < 2; ++k)
                    {
                        for (int l = 0; l < 2; ++l)
                        {
                            int index = (j - 1 + l * 2) % chromosomeSize;
                            while (index < 0)
                            {
                                index += chromosomeSize;
                            }

                            int next = matingPool.at(parents[k]).genes.at(index);

                            auto itr = find(edgeTable.at(matingPool.at(parents[k]).genes.at(j)).begin(), edgeTable.at(matingPool.at(parents[k]).genes.at(j)).end(), next);
                            if (itr == edgeTable.at(matingPool.at(parents[k]).genes.at(j)).end())
                            {
                                edgeTable.at(matingPool.at(parents[k]).genes.at(j)).push_back(next);
                            }
                        }
                    }
                }

                for (int j = 0; j < 2; ++j)
                {
                    Chromosome child;
                    auto table = edgeTable;

                    int current = matingPool.at(parents[j]).genes.at(0);
                    child.genes.push_back(current);

                    for (int k = 1; k < chromosomeSize; ++k)
                    {
                        for (auto& v : table)
                        {
                            auto itr = find(v.begin(), v.end(), current);
                            if (itr != v.end())
                            {
                                v.erase(itr);
                            }
                        }

                        if (table.at(current).size() > 0)
                        {
                            int min = table.at(current).at(0);
                            for (int l = 1; l < table.at(current).size(); ++l)
                            {
                                int m = table.at(table.at(current).at(l)).size();
                                if (m < table.at(min).size())
                                {
                                    min = table.at(current).at(l);
                                }
                            }

                            current = min;
                        }
                        else
                        {
                            for (int l = 0; l < chromosomeSize; ++l)
                            {
                                auto itr = find(child.genes.begin(), child.genes.end(), l);
                                if (itr == child.genes.end())
                                {
                                    current = l;
                                    break;
                                }
                            }
                        }

                        child.genes.push_back(current);
                    }

                    population.push_back(child);
                }
            }
            // 順序交叉
            else
            {
                int splitIndex = rand() % chromosomeSize;
                int splitLength = rand() % chromosomeSize;

                for (int j = 0; j < 2; ++j)
                {
                    Chromosome child;
                    child.genes = vector<int>(chromosomeSize, -1);

                    for (int k = 0; k < splitLength; ++k)
                    {
                        child.genes.at((splitIndex + k) % chromosomeSize) = matingPool.at(parents[j]).genes.at((splitIndex + k) % chromosomeSize);
                    }

                    int l = 0;
                    for (int m = 0; l < chromosomeSize - splitLength; ++m)
                    {
                        int n = matingPool.at(parents[1 - j]).genes.at((splitIndex + splitLength + m) % chromosomeSize);

                        auto itr = find(child.genes.begin(), child.genes.end(), n);
                        if (itr == child.genes.end())
                        {
                            child.genes.at((splitIndex + splitLength + l) % chromosomeSize) = n;
                            l++;
                        }
                    }

                    population.push_back(child);
                }
            }
        }
        else
        {
            population.push_back(matingPool.at(parents[0]));
            population.push_back(matingPool.at(parents[1]));
        }
    }

    // 個体数が奇数なら余った1個体はそのままコピーする
    if (populationSize % 2)
    {
        population.push_back(matingPool.at(populationSize - 1));
    }
}

// 突然変異
// 各個体に対しmutateProbabilityの確率で突然変異を行う
void GA::mutate()
{
    for (int i = 0; i < populationSize; ++i)
    {
        // ランダムに遺伝子を2つ選び交換する
        float p = (float)rand() / (float)RAND_MAX;

        if (p <= mutateProbability)
        {
            int a = rand() % chromosomeSize;
            int b = rand() % chromosomeSize;

            int t = population.at(i).genes.at(a);
            population.at(i).genes.at(a) = population.at(i).genes.at(b);
            population.at(i).genes.at(b) = t;
        }
    }
}

// アルゴリズムを実行する
void GA::solve()
{
    initialize();
    evaluate();

    for (int i = 1; i < generationSize - 1; ++i)
    {
        select();
        crossover();
        mutate();
        evaluate();
    }

    cout << "best solution: ";

    for (auto& g : best.genes)
    {
        cout << g << " ";
    }

    cout << "\nfitness: " << best.fitness << endl;
}

// ソート用の比較演算子
bool operator<(const GA::Chromosome &lhs, const GA::Chromosome &rhs)
{
    return lhs.fitness < rhs.fitness;
}