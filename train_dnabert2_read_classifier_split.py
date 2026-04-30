#!/usr/bin/env python3
import argparse
import os

import evaluate
import numpy as np
import pandas as pd
from datasets import Dataset
from transformers import (
    AutoModelForSequenceClassification,
    AutoTokenizer,
    DataCollatorWithPadding,
    Trainer,
    TrainingArguments,
)
from transformers.trainer_utils import get_last_checkpoint


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_csv", required=True)
    ap.add_argument("--test_csv", required=True)
    ap.add_argument("--model_name", default="zhihan1996/DNABERT-2-117M")
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--max_length", type=int, default=250)
    ap.add_argument("--min_class_count", type=int, default=20)
    ap.add_argument("--epochs", type=int, default=3)
    ap.add_argument("--batch_size", type=int, default=8)
    ap.add_argument("--eval_steps", type=int, default=1000)
    ap.add_argument("--save_steps", type=int, default=1000)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    train_df = pd.read_csv(args.train_csv)
    test_df = pd.read_csv(args.test_csv)

    train_df["label_taxid"] = train_df["label_taxid"].astype(str)
    test_df["label_taxid"] = test_df["label_taxid"].astype(str)

    counts = train_df["label_taxid"].value_counts()
    keep = set(counts[counts >= args.min_class_count].index)

    train_df = train_df[train_df["label_taxid"].isin(keep)].copy()
    test_df = test_df[test_df["label_taxid"].isin(keep)].copy()

    if len(train_df) == 0 or len(test_df) == 0:
        raise ValueError("No rows left after min_class_count filter")

    labels = sorted(train_df["label_taxid"].unique())
    label2id = {lab: i for i, lab in enumerate(labels)}
    id2label = {i: lab for lab, i in label2id.items()}

    train_df["label"] = train_df["label_taxid"].map(label2id)
    test_df["label"] = test_df["label_taxid"].map(label2id)

    print(f"Training rows after filter: {len(train_df)}")
    print(f"Test rows after filter: {len(test_df)}")
    print(f"Number of classes: {len(labels)}")

    tokenizer = AutoTokenizer.from_pretrained(args.model_name, trust_remote_code=True)
    model = AutoModelForSequenceClassification.from_pretrained(
        args.model_name,
        num_labels=len(labels),
        id2label=id2label,
        label2id=label2id,
        trust_remote_code=True,
    )

    train_ds = Dataset.from_pandas(train_df[["sequence", "label"]], preserve_index=False)
    test_ds = Dataset.from_pandas(test_df[["sequence", "label"]], preserve_index=False)

    def tokenize(batch):
        return tokenizer(
            batch["sequence"],
            truncation=True,
            padding=False,
            max_length=args.max_length,
        )

    train_ds = train_ds.map(tokenize, batched=True, num_proc=4)
    test_ds = test_ds.map(tokenize, batched=True, num_proc=4)

    acc = evaluate.load("accuracy")
    f1 = evaluate.load("f1")

    def compute_metrics(eval_pred):
        predictions, labels = eval_pred

        if isinstance(predictions, tuple):
            predictions = predictions[0]

        preds = np.argmax(predictions, axis=-1)

        return {
            "accuracy": acc.compute(predictions=preds, references=labels)["accuracy"],
            "macro_f1": f1.compute(
                predictions=preds,
                references=labels,
                average="macro",
            )["f1"],
        }

    training_args = TrainingArguments(
        output_dir=args.output_dir,
        evaluation_strategy="steps",
        eval_steps=args.eval_steps,
        save_strategy="steps",
        save_steps=args.save_steps,
        logging_strategy="steps",
        logging_steps=100,
        learning_rate=2e-5,
        per_device_train_batch_size=args.batch_size,
        per_device_eval_batch_size=args.batch_size,
        num_train_epochs=args.epochs,
        weight_decay=0.01,
        save_total_limit=2,
        load_best_model_at_end=True,
        metric_for_best_model="macro_f1",
        greater_is_better=True,
        fp16=False,
        dataloader_num_workers=4,
        report_to="none",
    )

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_ds,
        eval_dataset=test_ds,
        tokenizer=tokenizer,
        data_collator=DataCollatorWithPadding(tokenizer=tokenizer),
        compute_metrics=compute_metrics,
    )

    last_checkpoint = get_last_checkpoint(args.output_dir)
    if last_checkpoint is not None:
        print(f"Resuming from checkpoint: {last_checkpoint}")
    else:
        print("No checkpoint found; starting fresh.")

    trainer.train(resume_from_checkpoint=last_checkpoint)

    metrics = trainer.evaluate()
    print(metrics)

    trainer.save_model(args.output_dir)
    tokenizer.save_pretrained(args.output_dir)

    with open(os.path.join(args.output_dir, "label_map.tsv"), "w") as out:
        out.write("label_id\tlabel_taxid\n")
        for i, lab in id2label.items():
            out.write(f"{i}\t{lab}\n")

    with open(os.path.join(args.output_dir, "test_metrics.tsv"), "w") as out:
        for k, v in sorted(metrics.items()):
            out.write(f"{k}\t{v}\n")


if __name__ == "__main__":
    main()
