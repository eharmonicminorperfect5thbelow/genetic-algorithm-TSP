#ifndef INCLUDE_GA_HPP
#define INCLUDE_GA_HPP

#include <cmath>
#include <vector>

using namespace std;

class GA
{
public:
    // 都市
    // x座標, y座標
    struct City
    {
        double x;
        double y;

        City(double x, double y) : x(x), y(y) {}
    };

    // 染色体(chromosome)
    // 遺伝子(genes)とそれに対応する適合度(fitness)
    struct Chromosome
    {
        vector<int> genes;
        double fitness;
    };

    // コンストラクタ
    GA(const char* instance, int populationSize, int generationSize, float selectivePressure, float crossoverRatio, float crossoverProbability, float mutateProbability);
    
    // 計算を実行
    void solve();

private:
    // 都市間の距離を計算
    double dist(const City& a, const City& b)
    {
        return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
    }

    // インスタンスの都市の配列
    vector<City> cities;
    // 集団
    vector<Chromosome> population;
    // 選択(select)により得られた集団
    vector<Chromosome> matingPool;
    // その時点で最も適合度(fitness)の良い個体
    Chromosome best;
    // 最大の世代数
    int generationSize;
    // 集団内の個体数
    int populationSize;
    // 染色体の長さ
    int chromosomeSize;
    // 集団内での適合度の順位に対応した選択確率
    vector<float> selectProbabilities;
    // 交叉(crossover)でedge recombinationを用いる割合
    float crossoverRatio;
    // 交叉を行う確率
    float crossoverProbability;
    // 突然変異を行う確率
    float mutateProbability;

    // 初期集団の生成
    void initialize();
    // 適合度の計算
    void evaluate();
    // 選択
    void select();
    // 交叉
    void crossover();
    // 突然変異
    void mutate();
};

#endif
