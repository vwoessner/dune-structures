#ifndef DUNE_STRUCTURES_TIMECAPSULE_HH
#define DUNE_STRUCTURES_TIMECAPSULE_HH


template<typename F=double>
class TimeCapsule
{
  F time;
public:
  TimeCapsule ()
  {
    time = 0.0;
  }

  TimeCapsule (F t)
  {
    time = t;
  }

  // get time
  F getTime () const
  {
    return time;
  }

  // store time
  void setTime (F t)
  {
    time = t;
  }
};

#endif
