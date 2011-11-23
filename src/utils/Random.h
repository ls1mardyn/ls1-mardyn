class Random
{
 public:
   Random();
   void init(int seed);

   float rnd();

   int getIX() { return this->ix; }

 private:
   int ix, iy;
   float am;
};

