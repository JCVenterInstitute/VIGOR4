package org.jcvi.vigor.utils;

import org.junit.runners.Parameterized;
import org.junit.runners.model.RunnerScheduler;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

// From http://hwellmann.blogspot.com/2009/12/running-parameterized-junit-tests-in.html
public class ParallelizedParameterized extends Parameterized
{

    private static class ThreadPoolScheduler implements RunnerScheduler
    {
        private ExecutorService executor;

        public ThreadPoolScheduler()
        {
            int defaultThreadCount = Runtime.getRuntime().availableProcessors() * 2;
            String threads = System.getProperty("junit.parallel.threads");
            int numThreads = threads != null ? Integer.parseInt(threads): defaultThreadCount;
            executor = Executors.newFixedThreadPool(numThreads);
        }

        @Override
        public void finished()
        {
            executor.shutdown();
            try
            {
                executor.awaitTermination(5, TimeUnit.MINUTES);
            }
            catch (InterruptedException exc)
            {
                throw new RuntimeException(exc);
            }
        }

        @Override
        public void schedule(Runnable childStatement)
        {
            executor.submit(childStatement);
        }
    }

    public ParallelizedParameterized(Class klass) throws Throwable
    {
        super(klass);
        setScheduler(new ThreadPoolScheduler());
    }
}

