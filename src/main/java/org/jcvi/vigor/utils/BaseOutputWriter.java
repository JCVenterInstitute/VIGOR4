package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.exception.VigorException;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Consumer;

public abstract class BaseOutputWriter implements IOutputWriter, IConfigurable {

    private static Logger LOGGER = LogManager.getLogger(BaseOutputWriter.class);
    private static Consumer<BufferedWriter> nullHandler = (b) -> {};

    protected boolean multiFile = false;

    public void configure(VigorConfiguration config) {
        multiFile = config.getOrDefault(ConfigurationParameters.MultiFile, false);
    }

    public Consumer<BufferedWriter> getOnOpen() {
        return nullHandler;
    }

    public Consumer<BufferedWriter> getOnClose() {
        return nullHandler;
    }

    public WriterBundle getWriter(Outfiles outfiles,
                                  OutputContext context,
                                  OutputContext.Key ... toClose )
            throws IOException, VigorException
    {
        List<BufferedWriter> writers = new ArrayList<>();
        List<BufferedWriter> borrowedBuffers = new ArrayList<>();
        Path baseFilePath = outfiles.getBaseFilePath(getExtension());
        writers.add(outfiles.getWriter(baseFilePath, getOnOpen(), getOnClose()));
        borrowedBuffers.addAll(writers);

        if (multiFile) {
            EnumSet<OutputContext.Key> closeSet = toClose.length == 0 ? EnumSet.noneOf(OutputContext.Key.class): EnumSet.copyOf(Arrays.asList(toClose));
            BufferedWriter bw;
            for (OutputContext.Key key : context.keySet()) {
                if (!context.getContext(key).isPresent()) {
                    LOGGER.warn("No context value provided for key {}", key);
                    continue;
                }
                String value = context.getContext(key).get();
                switch (key) {
                    case GENOME:
                        String sequenceFile = OutputWriterUtils.getSequenceFilePath(value) + "." +  getExtension();
                        bw = outfiles.getWriter(Paths.get(sequenceFile), getOnOpen(), getOnClose());
                        break;
                    case GENE:
                    case PEP:
                        String geneFilePath = OutputWriterUtils.getGeneFilePath(value) + "." + getExtension();
                        bw = outfiles.getWriter(Paths.get(geneFilePath), getOnOpen(), getOnClose());
                        break;
                    default:
                        throw new VigorException("Unexpected context key " + key.toString());
                }
                writers.add(bw);
                if (!closeSet.contains(key)) {
                    borrowedBuffers.add(bw);
                }
            }
        }
        WriterBundle bundle = new WriterBundle(writers.toArray(new BufferedWriter[]{}));
        for (BufferedWriter writer: borrowedBuffers) {
            bundle.borrowed(writer);
        }
        return bundle;
    }

}
