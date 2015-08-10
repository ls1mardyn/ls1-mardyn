class Random
{
 public:
   Random();
   void init(int seed);

   float rnd();

   int getIX() { return this->ix_muVT; }

 private:
   int ix_muVT, iy_muVT;
   float am_muVT;
};

